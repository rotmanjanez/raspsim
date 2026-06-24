#include "x86sim-support/syscall-linux.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <expected>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using x86sim::address_t;
using x86sim::word_t;

constexpr address_t stack_top = 0x7fe000000000ULL;
constexpr std::uint64_t stack_size = 1024 * 1024;

constexpr unsigned char elf_class_64 = 2;
constexpr unsigned char elf_data_lsb = 1;
constexpr std::uint16_t et_exec = 2;
constexpr std::uint16_t em_x86_64 = 62;
constexpr std::uint32_t pt_load = 1;
constexpr std::uint32_t pt_interp = 3;
constexpr std::uint32_t pf_x = 1;
constexpr std::uint32_t pf_w = 2;
constexpr std::uint32_t pf_r = 4;

struct Elf64Ehdr {
  std::array<unsigned char, 16> e_ident;
  std::uint16_t e_type;
  std::uint16_t e_machine;
  std::uint32_t e_version;
  std::uint64_t e_entry;
  std::uint64_t e_phoff;
  std::uint64_t e_shoff;
  std::uint32_t e_flags;
  std::uint16_t e_ehsize;
  std::uint16_t e_phentsize;
  std::uint16_t e_phnum;
  std::uint16_t e_shentsize;
  std::uint16_t e_shnum;
  std::uint16_t e_shstrndx;
};

struct Elf64Phdr {
  std::uint32_t p_type;
  std::uint32_t p_flags;
  std::uint64_t p_offset;
  std::uint64_t p_vaddr;
  std::uint64_t p_paddr;
  std::uint64_t p_filesz;
  std::uint64_t p_memsz;
  std::uint64_t p_align;
};

static_assert(sizeof(Elf64Ehdr) == 64);
static_assert(sizeof(Elf64Phdr) == 56);

struct LoadedElf {
  address_t entry = 0;
  address_t phdr = 0;
  word_t phent = 0;
  word_t phnum = 0;
};

[[nodiscard]] bool has_slash(std::string_view value) {
  return value.find('/') != std::string_view::npos;
}

[[nodiscard]] bool executable_candidate(const std::filesystem::path& path) {
  std::error_code ec;
  return std::filesystem::is_regular_file(path, ec);
}

[[nodiscard]] std::filesystem::path resolve_program(std::string_view program) {
  std::filesystem::path requested(program);
  if (has_slash(program)) {
    if (executable_candidate(requested))
      return requested;
    throw std::runtime_error("program not found: " + std::string(program));
  }

  if (executable_candidate(requested))
    return requested;

  const char* path_env = std::getenv("PATH");
  if (!path_env)
    throw std::runtime_error("program not found: " + std::string(program));

  std::string_view paths(path_env);
  while (true) {
    const auto separator = paths.find(':');
    std::string_view entry = paths.substr(0, separator);
    std::filesystem::path directory = entry.empty() ? "." : std::filesystem::path(entry);
    std::filesystem::path candidate = directory / requested;
    if (executable_candidate(candidate))
      return candidate;

    if (separator == std::string_view::npos)
      break;
    paths.remove_prefix(separator + 1);
  }

  throw std::runtime_error("program not found in PATH: " + std::string(program));
}

[[nodiscard]] std::vector<std::byte> read_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input)
    throw std::runtime_error("cannot open " + path.string());

  const auto size = input.tellg();
  if (size < 0)
    throw std::runtime_error("cannot determine size of " + path.string());

  std::vector<std::byte> bytes(static_cast<std::size_t>(size));
  input.seekg(0);
  if (!bytes.empty())
    input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  if (!input)
    throw std::runtime_error("cannot read " + path.string());
  return bytes;
}

template<typename T>
[[nodiscard]] const T& read_struct(std::span<const std::byte> bytes, std::uint64_t offset) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(T))
    throw std::runtime_error("truncated ELF file");

  const auto* ptr = bytes.data() + offset;
  return *reinterpret_cast<const T*>(ptr);
}

[[nodiscard]] address_t page_floor(address_t value) {
  return value & ~(x86sim::Machine::kPageSize - 1);
}

[[nodiscard]] address_t page_ceil(address_t value) {
  const address_t mask = x86sim::Machine::kPageSize - 1;
  if (value > std::numeric_limits<address_t>::max() - mask)
    throw std::runtime_error("address range overflow");
  return (value + mask) & ~mask;
}

[[nodiscard]] x86sim::Protection protection_from_elf(std::uint32_t flags) {
  x86sim::Protection protection = x86sim::Protection::none;
  if ((flags & pf_r) != 0)
    protection = protection | x86sim::Protection::read;
  if ((flags & pf_w) != 0)
    protection = protection | x86sim::Protection::write;
  if ((flags & pf_x) != 0)
    protection = protection | x86sim::Protection::execute;
  return protection;
}

void checked_map(x86sim::Machine& machine, address_t start, std::uint64_t size, x86sim::Protection protection) {
  auto mapped = machine.map(start, size, protection);
  if (!mapped)
    throw std::runtime_error("cannot map guest memory");
}

void checked_write(x86sim::Machine& machine, address_t start, std::span<const std::byte> bytes) {
  auto written = machine.write_memory(start, bytes);
  if (!written)
    throw std::runtime_error("cannot write guest memory");
}

[[nodiscard]] LoadedElf load_elf(x86sim::Machine& machine, const std::filesystem::path& path) {
  const auto file = read_file(path);
  const Elf64Ehdr& header = read_struct<Elf64Ehdr>(file, 0);

  if (header.e_ident[0] != 0x7f || header.e_ident[1] != 'E' || header.e_ident[2] != 'L' || header.e_ident[3] != 'F')
    throw std::runtime_error(path.string() + " is not an ELF file");
  if (header.e_ident[4] != elf_class_64 || header.e_ident[5] != elf_data_lsb)
    throw std::runtime_error(path.string() + " is not a little-endian ELF64 file");
  if (header.e_type != et_exec)
    throw std::runtime_error(path.string() + " is not an ET_EXEC executable; PIE/dynamic loading is not supported yet");
  if (header.e_machine != em_x86_64)
    throw std::runtime_error(path.string() + " is not an x86-64 ELF executable");
  if (header.e_phentsize != sizeof(Elf64Phdr))
    throw std::runtime_error(path.string() + " has unsupported program-header size");

  if (header.e_phoff > file.size() || header.e_phnum > (file.size() - header.e_phoff) / sizeof(Elf64Phdr))
    throw std::runtime_error(path.string() + " has truncated program headers");

  for (std::uint16_t i = 0; i < header.e_phnum; ++i) {
    const Elf64Phdr& phdr = read_struct<Elf64Phdr>(file, header.e_phoff + static_cast<std::uint64_t>(i) * sizeof(Elf64Phdr));
    if (phdr.p_type == pt_interp)
      throw std::runtime_error(path.string() + " uses PT_INTERP; dynamic loading is not supported yet");
  }

  LoadedElf loaded{
      .entry = header.e_entry,
      .phdr = 0,
      .phent = header.e_phentsize,
      .phnum = header.e_phnum,
  };

  for (std::uint16_t i = 0; i < header.e_phnum; ++i) {
    const Elf64Phdr& phdr = read_struct<Elf64Phdr>(file, header.e_phoff + static_cast<std::uint64_t>(i) * sizeof(Elf64Phdr));
    if (phdr.p_type != pt_load || phdr.p_memsz == 0)
      continue;
    if (phdr.p_filesz > phdr.p_memsz)
      throw std::runtime_error(path.string() + " has invalid segment sizes");
    if (phdr.p_offset > file.size() || phdr.p_filesz > file.size() - phdr.p_offset)
      throw std::runtime_error(path.string() + " has a truncated load segment");
    if (phdr.p_vaddr > std::numeric_limits<address_t>::max() - phdr.p_memsz)
      throw std::runtime_error(path.string() + " has an overflowing load segment");

    const address_t map_start = page_floor(phdr.p_vaddr);
    const address_t map_end = page_ceil(phdr.p_vaddr + phdr.p_memsz);
    checked_map(machine, map_start, map_end - map_start, protection_from_elf(phdr.p_flags));

    if (phdr.p_filesz != 0) {
      checked_write(machine, phdr.p_vaddr,
                    std::span<const std::byte>(file.data() + phdr.p_offset, static_cast<std::size_t>(phdr.p_filesz)));
    }

    if (header.e_phoff >= phdr.p_offset && header.e_phoff < phdr.p_offset + phdr.p_filesz)
      loaded.phdr = phdr.p_vaddr + (header.e_phoff - phdr.p_offset);
  }

  return loaded;
}

void append_word(std::vector<std::byte>& bytes, word_t value) {
  for (unsigned i = 0; i < 8; ++i)
    bytes.push_back(static_cast<std::byte>((value >> (i * 8)) & 0xff));
}

[[nodiscard]] address_t setup_stack(x86sim::Machine& machine, std::span<const std::string> argv, const LoadedElf& loaded) {
  const address_t stack_base = stack_top - stack_size;
  checked_map(machine, stack_base, stack_size, x86sim::Protection::read | x86sim::Protection::write);

  address_t cursor = stack_top;
  std::vector<address_t> argv_addresses(argv.size());
  for (std::size_t i = argv.size(); i > 0; --i) {
    const std::string& arg = argv[i - 1];
    if (arg.size() + 1 > cursor - stack_base)
      throw std::runtime_error("argument stack overflow");

    cursor -= arg.size() + 1;
    std::vector<std::byte> bytes(arg.size() + 1);
    std::memcpy(bytes.data(), arg.data(), arg.size());
    checked_write(machine, cursor, bytes);
    argv_addresses[i - 1] = cursor;
  }

  cursor -= 16;
  constexpr std::array<std::byte, 16> random_bytes{
      std::byte{0x52}, std::byte{0x41}, std::byte{0x53}, std::byte{0x50}, std::byte{0x73}, std::byte{0x69}, std::byte{0x6d},
      std::byte{0x21}, std::byte{0x10}, std::byte{0x32}, std::byte{0x54}, std::byte{0x76}, std::byte{0x98}, std::byte{0xba},
      std::byte{0xdc}, std::byte{0xfe},
  };
  checked_write(machine, cursor, random_bytes);
  const address_t random_address = cursor;

  cursor &= ~static_cast<address_t>(0xf);

  constexpr word_t at_null = 0;
  constexpr word_t at_phdr = 3;
  constexpr word_t at_phent = 4;
  constexpr word_t at_phnum = 5;
  constexpr word_t at_pagesz = 6;
  constexpr word_t at_base = 7;
  constexpr word_t at_flags = 8;
  constexpr word_t at_entry = 9;
  constexpr word_t at_uid = 11;
  constexpr word_t at_euid = 12;
  constexpr word_t at_gid = 13;
  constexpr word_t at_egid = 14;
  constexpr word_t at_clktck = 17;
  constexpr word_t at_secure = 23;
  constexpr word_t at_random = 25;

  std::vector<word_t> words;
  words.push_back(argv.size());
  for (address_t arg_address : argv_addresses)
    words.push_back(arg_address);
  words.push_back(0); // argv terminator
  words.push_back(0); // envp terminator
  if (loaded.phdr != 0) {
    words.push_back(at_phdr);
    words.push_back(loaded.phdr);
    words.push_back(at_phent);
    words.push_back(loaded.phent);
    words.push_back(at_phnum);
    words.push_back(loaded.phnum);
  }
  words.push_back(at_pagesz);
  words.push_back(x86sim::Machine::kPageSize);
  words.push_back(at_base);
  words.push_back(0);
  words.push_back(at_flags);
  words.push_back(0);
  words.push_back(at_entry);
  words.push_back(loaded.entry);
  words.push_back(at_uid);
  words.push_back(1000);
  words.push_back(at_euid);
  words.push_back(1000);
  words.push_back(at_gid);
  words.push_back(1000);
  words.push_back(at_egid);
  words.push_back(1000);
  words.push_back(at_clktck);
  words.push_back(100);
  words.push_back(at_secure);
  words.push_back(0);
  words.push_back(at_random);
  words.push_back(random_address);
  words.push_back(at_null);
  words.push_back(0);

  if (words.size() * sizeof(word_t) > cursor - stack_base)
    throw std::runtime_error("argument stack overflow");
  cursor -= words.size() * sizeof(word_t);
  cursor &= ~static_cast<address_t>(0xf);

  std::vector<std::byte> stack_words;
  stack_words.reserve(words.size() * sizeof(word_t));
  for (word_t word : words)
    append_word(stack_words, word);
  checked_write(machine, cursor, stack_words);

  return cursor;
}

void usage(const char* argv0) {
  std::cerr << "usage: " << argv0 << " [--core=seq|ooo] [--seq|--ooo] [--] PROGRAM [ARG...]\n";
}

[[nodiscard]] x86sim::CoreModel parse_core(std::string_view value) {
  if (value == "seq" || value == "sequential")
    return x86sim::CoreModel::sequential;
  if (value == "ooo" || value == "out-of-order" || value == "out_of_order")
    return x86sim::CoreModel::out_of_order;
  throw std::runtime_error("unknown core model: " + std::string(value));
}

} // namespace

int main(int argc, char** argv) {
  x86sim::CoreModel core = x86sim::CoreModel::sequential;
  int program_index = 1;
  for (; program_index < argc; ++program_index) {
    const std::string_view arg(argv[program_index]);
    if (arg == "--") {
      ++program_index;
      break;
    }
    if (arg == "-h" || arg == "--help") {
      usage(argv[0]);
      return 0;
    }
    if (arg == "--seq") {
      core = x86sim::CoreModel::sequential;
      continue;
    }
    if (arg == "--ooo") {
      core = x86sim::CoreModel::out_of_order;
      continue;
    }
    constexpr std::string_view core_prefix = "--core=";
    if (arg.starts_with(core_prefix)) {
      core = parse_core(arg.substr(core_prefix.size()));
      continue;
    }
    if (arg.starts_with("--")) {
      std::cerr << "x86sim-linux: unknown option: " << arg << "\n";
      usage(argv[0]);
      return 2;
    }
    break;
  }

  if (program_index >= argc) {
    usage(argv[0]);
    return 2;
  }

  try {
    const std::filesystem::path program_path = resolve_program(argv[program_index]);

    int exit_status = 0;
    auto callbacks = x86sim::linux_syscalls::host(
        x86sim::linux_syscalls::SysRead{} | x86sim::linux_syscalls::SysReadlink{} |
        x86sim::linux_syscalls::SysWrite{} | x86sim::linux_syscalls::SysPreadPwrite{} |
        x86sim::linux_syscalls::SysOpen{} |
        x86sim::linux_syscalls::SysClose{} | x86sim::linux_syscalls::SysFstat{} |
        x86sim::linux_syscalls::SysStat{} |
        x86sim::linux_syscalls::SysFileSystem{} |
        x86sim::linux_syscalls::SysLseek{} |
        x86sim::linux_syscalls::SysIoctl{} | x86sim::linux_syscalls::SysFcntl{} |
        x86sim::linux_syscalls::SysSocket{} |
        x86sim::linux_syscalls::SysConnect{} |
        x86sim::linux_syscalls::SysSignals{} |
        x86sim::linux_syscalls::SysBrk{} |
        x86sim::linux_syscalls::SysArchPrctl{} | x86sim::linux_syscalls::SysGetIdentity{} |
        x86sim::linux_syscalls::SysUname{} | x86sim::linux_syscalls::SysGetcwd{} |
        x86sim::linux_syscalls::SysSetTidAddress{} |
        x86sim::linux_syscalls::SysSetRobustList{} | x86sim::linux_syscalls::SysRseq{} |
        x86sim::linux_syscalls::SysPrlimit64{} |
        x86sim::linux_syscalls::SysMmap{} | x86sim::linux_syscalls::SysMunmap{} |
        x86sim::linux_syscalls::SysExit{&exit_status} | x86sim::linux_syscalls::SysExitGroup{&exit_status});

    x86sim::Options options;
    options.core = core;
    options.quiet = true;
    options.log.log_filename.clear();
    options.log.loglevel = 100;
    options.log.log_on_console = false;

    x86sim::Machine machine(callbacks, options);
    std::vector<std::string> guest_argv;
    guest_argv.reserve(static_cast<std::size_t>(argc - program_index));
    guest_argv.emplace_back(argv[program_index]);
    for (int i = program_index + 1; i < argc; ++i)
      guest_argv.emplace_back(argv[i]);

    const LoadedElf loaded = load_elf(machine, program_path);
    const address_t rsp = setup_stack(machine, guest_argv, loaded);

    x86sim::RegisterFile& registers = machine.register_file(0);
    registers[x86sim::Register::rip] = loaded.entry;
    registers[x86sim::Register::rsp] = rsp;

    const x86sim::RunResult result = machine.run();
    switch (result.reason) {
    case x86sim::StopReason::guest_exit:
      return exit_status;
    case x86sim::StopReason::unsupported_syscall:
      std::cerr << "x86sim-linux: unsupported syscall";
      if (!result.message.empty())
        std::cerr << ": " << result.message;
      std::cerr << "\n";
      return 1;
    case x86sim::StopReason::x86_exception:
      std::cerr << "x86sim-linux: x86 exception: " << result.message << "\n";
      return 1;
    case x86sim::StopReason::instruction_limit:
      std::cerr << "x86sim-linux: instruction limit reached\n";
      return 1;
    case x86sim::StopReason::host_request:
      std::cerr << "x86sim-linux: stopped by host request";
      if (!result.message.empty())
        std::cerr << ": " << result.message;
      std::cerr << "\n";
      return 1;
    }
  } catch (const std::exception& error) {
    std::cerr << "x86sim-linux: " << error.what() << "\n";
    return 1;
  }

  return 1;
}
