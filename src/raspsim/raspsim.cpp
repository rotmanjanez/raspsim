//
// PTLsim: Cycle Accurate x86-64 Simulator
// RASPsim application
//
// Copyright 2020-2020 Alexis Engelke <engelke@in.tum.de>
//
#include "x86sim-support/cpuid.hpp"
#include "x86sim/addrspace.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"
#include "x86sim/logging.hpp"

#include "config.h"

#include <cerrno>
#include <charconv>
#include <cstddef>
#include <fstream>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

extern ConfigurationParser<x86sim::Options> configparser;

namespace {

namespace logging = x86sim::logging;

using std::operator""sv;

void print_hex_bytes(std::span<const std::byte> bytes, size_t splitat = 16) {
  for (size_t i = 0; i < bytes.size(); i++) {
    logging::eprint("{:02x}", std::to_integer<unsigned>(bytes[i]));
    if (((i % splitat) == (splitat - 1)) && (i != bytes.size() - 1))
      logging::eprint("\n");
    else if (i != bytes.size() - 1)
      logging::eprint(" ");
  }
  logging::eprintln("");
}

[[nodiscard]] std::optional<x86sim::word_t> parse_hex(std::string_view str) {
  x86sim::word_t value;
  if (str.starts_with("0x"))
    str = str.substr(2);
  auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), value, 16);
  if (ec == std::errc())
    return value;
  return std::nullopt;
}

[[nodiscard]] std::optional<x86sim::Register> register_from_name(std::string_view name) {
  using enum x86sim::Register;
  static constexpr std::pair<std::string_view, x86sim::Register> names[] = {
      {"rax", rax}, {"rcx", rcx}, {"rdx", rdx}, {"rbx", rbx}, {"rsp", rsp}, {"rbp", rbp},
      {"rsi", rsi}, {"rdi", rdi}, {"r8", r8},   {"r9", r9},   {"r10", r10}, {"r11", r11},
      {"r12", r12}, {"r13", r13}, {"r14", r14}, {"r15", r15}, {"rip", rip}, {"flags", flags},
  };

  for (auto [reg_name, reg] : names) {
    if (name == reg_name)
      return reg;
  }
  return std::nullopt;
}

[[nodiscard]] std::optional<std::pair<x86sim::XmmRegister, bool>> xmm_half_from_name(std::string_view name) {
  bool high = false;
  if (name.starts_with("xmml"))
    name = name.substr(4);
  else if (name.starts_with("xmmh")) {
    high = true;
    name = name.substr(4);
  } else {
    return std::nullopt;
  }

  unsigned value = 0;
  auto [ptr, ec] = std::from_chars(name.data(), name.data() + name.size(), value, 10);
  if (ec != std::errc() || ptr != name.data() + name.size() || value > 15)
    return std::nullopt;

  const auto reg = static_cast<x86sim::XmmRegister>(static_cast<int>(x86sim::XmmRegister::xmm0) + value * 2);
  return std::pair{reg, high};
}

[[nodiscard]] std::optional<x86sim::Protection> protection_from_string(std::string_view prot) {
  using enum x86sim::Protection;
  if (prot == "ro")
    return read;
  if (prot == "rw")
    return read | write;
  if (prot == "rx")
    return read | execute;
  if (prot == "rwx")
    return read | write | execute;
  return std::nullopt;
}

// raspsim is a register/memory-level simulator: a 0x80 interrupt signals the
// guest is done, and any other syscall just returns -ENOSYS so execution can
// continue (or be inspected).
class CliHost : public x86sim::HostCallbacks {
public:
  [[nodiscard]] x86sim::CpuidResult cpuid(x86sim::Machine&, x86sim::CpuState&, x86sim::AddressSpace&,
                                          x86sim::CpuidRequest request) noexcept override {
    return x86sim::defaults::default_cpuid(request);
  }

  [[nodiscard]] x86sim::SyscallResult syscall(x86sim::Machine&, x86sim::CpuState& context, x86sim::AddressSpace&,
                                              x86sim::SyscallKind kind) override {
    using x86sim::Register;
    using x86sim::StopReason;
    if (kind == x86sim::SyscallKind::int80)
      return {.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};

    context[Register::rax] = static_cast<x86sim::word_t>(-ENOSYS);
    if (kind == x86sim::SyscallKind::syscall64)
      context[Register::rip] = context[Register::rcx];

    return {.reason = StopReason::host_request, .continue_execution = true, .message = {}};
  }
};

[[nodiscard]] std::optional<std::vector<std::byte>> parse_hex_bytes(std::string_view hex) {
  if ((hex.size() & 1) != 0)
    return std::nullopt;

  std::vector<std::byte> bytes(hex.size() / 2);
  for (size_t i = 0; i < bytes.size(); ++i) {
    unsigned value = 0;
    auto byte = hex.substr(i * 2, 2);
    auto [ptr, ec] = std::from_chars(byte.data(), byte.data() + byte.size(), value, 16);
    if (ec != std::errc())
      return std::nullopt;
    bytes[i] = static_cast<std::byte>(value);
  }

  return bytes;
}

bool handle_config_arg(x86sim::AddressSpace* address_space, x86sim::RegisterFile* registers, x86sim::Options& options,
                       std::string_view line, std::vector<x86sim::address_t>* dump_pages) {
  if (line.empty())
    return false;

  std::vector<std::string_view> toks;
  for (auto tok : std::views::split(line, " "sv)) {
    auto sv = std::string_view(tok.begin(), tok.end());
    if (!sv.empty())
      toks.push_back(sv);
  }
  if (toks.empty() || toks[0][0] == '#')
    return false;

  if (toks[0] == "Fnox87") {
    options.x87 = false;
  } else if (toks[0] == "Fnosse") {
    options.sse = false;
  } else if (toks[0] == "Fnocache") {
    options.debug.perfect_cache = true;
  } else if (toks[0] == "Fstbrpred") {
    options.debug.static_branchpred = true;
  } else if (address_space == nullptr || registers == nullptr) {
    return false;
  } else if (toks[0][0] == 'M') {
    if (toks.size() != 2) {
      logging::eprintln("Error: option {} has wrong number of arguments", line);
      return true;
    }
    auto addr = parse_hex(toks[0].substr(1));
    auto prot = protection_from_string(toks[1]);
    if (!addr || (*addr & (x86sim::AddressSpace::kPageSize - 1)) != 0) {
      logging::eprintln("Error: invalid value {}", toks[0]);
      return true;
    }
    if (!prot) {
      logging::eprintln("Error: invalid mem prot {}", toks[1]);
      return true;
    }
    if (auto result = address_space->map(*addr, x86sim::AddressSpace::kPageSize, *prot); !result) {
      logging::eprintln("Error: {}", result.error());
      return true;
    }
  } else if (toks[0][0] == 'W') {
    if (toks.size() != 2) {
      logging::eprintln("Error: option {} has wrong number of arguments", line);
      return true;
    }
    auto addr = parse_hex(toks[0].substr(1));
    auto bytes = parse_hex_bytes(toks[1]);
    if (!addr) {
      logging::eprintln("Error: invalid value {}", toks[0]);
      return true;
    }
    if (!bytes || bytes->size() > x86sim::AddressSpace::kPageSize - (*addr & (x86sim::AddressSpace::kPageSize - 1))) {
      logging::eprintln("Error: arg has odd size or crosses page boundary {}", (void*)*addr);
      return true;
    }
    if (auto result = address_space->write(*addr, std::as_bytes(std::span(*bytes))); !result) {
      logging::eprintln("Error: {}", result.error());
      return true;
    }
  } else if (toks[0][0] == 'D') {
    if (toks.size() != 1) {
      logging::eprintln("Error: option {} has wrong number of arguments", line);
      return true;
    }
    auto addr = parse_hex(toks[0].substr(1));
    if (!addr) {
      logging::eprintln("Error: invalid value {}", toks[0]);
      return true;
    }
    dump_pages->push_back(*addr & ~(x86sim::AddressSpace::kPageSize - 1));
  } else {
    if (toks.size() != 2) {
      logging::eprintln("Error: option {} has wrong number of arguments", line);
      return true;
    }
    auto value = parse_hex(toks[1]);
    if (!value) {
      logging::eprintln("Error: invalid value {}", toks[1]);
      return true;
    }
    if (auto reg = register_from_name(toks[0])) {
      (*registers)[*reg] = *value;
    } else if (auto xmm = xmm_half_from_name(toks[0])) {
      auto current = static_cast<x86sim::XmmValue>((*registers)[xmm->first]);
      if (xmm->second)
        current.hi = *value;
      else
        current.lo = *value;
      (*registers)[xmm->first] = current;
    } else {
      logging::eprintln("Error: invalid register {}", toks[0]);
      return true;
    }
  }

  return false;
}

} // namespace

int main(int argc, char** argv) {
  configparser.setup();

  x86sim::Options options;
  int ptlsim_arg_count = 1 + configparser.parse(options, argc - 1, argv + 1);
  if (ptlsim_arg_count == 0)
    ptlsim_arg_count = argc;

  std::vector<std::string> commands;
  for (unsigned i = ptlsim_arg_count; i < static_cast<unsigned>(argc); i++) {
    if (argv[i][0] == '@') {
      std::ifstream is(argv[i] + 1);
      if (!is) {
        logging::eprintln("Warning: cannot open command list file '{}'", argv[i]);
        continue;
      }
      std::string line;
      while (std::getline(is, line)) {
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
          line.erase(comment_pos);
        commands.push_back(std::move(line));
      }
    } else {
      commands.emplace_back(argv[i]);
    }
  }

  for (const std::string& command : commands)
    handle_config_arg(nullptr, nullptr, options, command, nullptr);

  CliHost host;
  x86sim::Machine machine(host, options);
  x86sim::AddressSpace address_space(machine);
  x86sim::CpuState initial_state;
  x86sim::RegisterFile& registers = initial_state;
  std::vector<x86sim::address_t> dump_pages;

  bool parse_err = false;
  for (const std::string& command : commands)
    parse_err |= handle_config_arg(&address_space, &registers, options, command, &dump_pages);

  if (parse_err) {
    logging::eprintln("Error: could not parse all arguments");
    return 1;
  }

  logging::eprintln("\n=== Switching to simulation mode at rip {} ===\n", registers[x86sim::Register::rip]);
  logging::eprintln("Baseline state:\n{}", registers);

  x86sim::RunResult result = machine.run(initial_state, address_space);
  x86sim::RegisterFile& final_registers = initial_state;

  if (result.reason == x86sim::StopReason::x86_exception) {
    logging::eprintln("{}", result);
    return 1;
  }

  logging::eprintln("End state:\n{}", final_registers);

  for (x86sim::address_t addr : dump_pages) {
    std::vector<std::byte> page(x86sim::AddressSpace::kPageSize);
    if (auto read = address_space.read(addr, page)) {
      logging::eprintln("Dump of memory at {}:", (void*)addr);
      print_hex_bytes(page);
    } else {
      logging::eprintln("Error dumping memory at {}: {}", (void*)addr, read.error());
    }
  }

  logging::eprintln("\n=== Exiting after full simulation at rip {} ({}) ===\n", final_registers[x86sim::Register::rip],
                    result.stats);

  return 0;
}
