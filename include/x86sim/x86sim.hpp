#ifndef RASPSIM_RASPSIM_HPP
#define RASPSIM_RASPSIM_HPP

#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "x86sim/registerfile.hpp"

namespace x86sim {

class AddressSpace;
struct Context;
struct MachineImpl;

using address_t = std::uint64_t;
enum class Protection : std::uint8_t { none = 0, read = 1, write = 2, execute = 4 };
enum class CoreModel { out_of_order, sequential };

enum class SyscallKind { int80, syscall64, sysenter };
enum class StopReason { guest_exit, instruction_limit, x86_exception, host_request, unsupported_syscall };
enum class MemoryError { unaligned_address, zero_size, unmapped_address, out_of_memory, mapping_failed };

[[nodiscard]] constexpr Protection operator|(Protection lhs, Protection rhs) noexcept {
  return static_cast<Protection>(static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}

[[nodiscard]] constexpr Protection operator&(Protection lhs, Protection rhs) noexcept {
  return static_cast<Protection>(static_cast<std::uint8_t>(lhs) & static_cast<std::uint8_t>(rhs));
}

[[nodiscard]] constexpr bool has_protection(Protection value, Protection flag) noexcept {
  return (static_cast<std::uint8_t>(value) & static_cast<std::uint8_t>(flag)) == static_cast<std::uint8_t>(flag);
}


struct CpuidRequest {
  std::uint32_t function = 0;
  std::uint32_t subfunction = 0;
};

struct CpuidResult {
  std::uint32_t eax = 0;
  std::uint32_t ebx = 0;
  std::uint32_t ecx = 0;
  std::uint32_t edx = 0;
};

struct SyscallResult {
  StopReason reason = StopReason::host_request;
  bool continue_execution = false;
  std::string message;
};

class Machine;

class HostCallbacks {
public:
  virtual ~HostCallbacks() = default;
  virtual SyscallResult syscall(Machine&, Context&, SyscallKind) noexcept = 0;
  virtual CpuidResult cpuid(Machine&, Context&, CpuidRequest) noexcept = 0;
};

struct Options {
  static constexpr address_t invalid_rip = std::numeric_limits<address_t>::max();

  CoreModel core = CoreModel::out_of_order;
  std::uint64_t core_count = 1;
  bool sse = true;
  bool x87 = true;

  bool quiet = false;
  std::filesystem::path log_filename = "ptlsim.log";
  std::uint64_t loglevel = 0;
  std::uint64_t start_log_at_iteration = std::numeric_limits<std::uint64_t>::max();
  address_t start_log_at_rip = invalid_rip;
  bool log_on_console = false;
  bool log_ptlsim_boot = false;
  std::uint64_t log_buffer_size = 524288;
  std::string mm_logfile;
  std::uint64_t mm_log_buffer_size = 16384;
  bool enable_inline_mm_logging = false;
  bool enable_mm_validate = false;

  bool event_log_enabled = false;
  std::uint64_t event_log_ring_buffer_size = 32768;
  bool flush_event_log_every_cycle = false;
  address_t log_backwards_from_trigger_rip = invalid_rip;
  bool dump_state_now = false;
  bool abort_at_end = false;

  address_t start_at_rip = 0;
  bool include_dyn_linker = false;
  bool trigger_mode = false;
  std::uint64_t pause_at_startup = 0;

  std::uint64_t stop_at_user_insns = std::numeric_limits<std::uint64_t>::max();
  std::uint64_t stop_at_iteration = std::numeric_limits<std::uint64_t>::max();
  address_t stop_at_rip = invalid_rip;

  bool perfect_cache = false;
  bool static_branchpred = false;

  std::string bbcache_dump_filename;
  std::uint64_t sequential_mode_insns = 0;
  bool exit_after_fullsim = false;
};

struct RunOptions {
  std::optional<std::uint64_t> instruction_limit = std::nullopt;
};

struct Stats {
  std::uint64_t cycles = 0;
  std::uint64_t instructions = 0;
};

struct X86Exception {
  std::uint8_t vector = 0;
  std::uint32_t error_code = 0;
  address_t virtual_address = 0;
  address_t rip = 0;
  std::string message;
  std::string context;
};

struct RunResult {
  StopReason reason;
  Stats stats;
  std::optional<X86Exception> x86_exception;
  std::string message;
};


class Machine {
public:
  static constexpr std::uint64_t kPageSize = 4096;

  explicit Machine(HostCallbacks&, Options = {});
  ~Machine();

  [[nodiscard]] std::expected<void, MemoryError> map(address_t start, std::uint64_t size, Protection) noexcept;
  void unmap(address_t start, std::uint64_t size) noexcept;
  [[nodiscard]] std::expected<std::span<const std::byte>, MemoryError> read_memory(address_t start,
                                                                                  std::size_t num_bytes) const noexcept;
  [[nodiscard]] std::expected<void, MemoryError> write_memory(address_t start, std::span<const std::byte>) noexcept;

  [[nodiscard]] RunResult run(RunOptions = {});

  [[nodiscard]] const Options& options() const noexcept;
  [[nodiscard]] Stats stats() const noexcept;
  [[nodiscard]] RegisterFile& register_file(std::size_t core_index = 0) noexcept;
  [[nodiscard]] const RegisterFile& register_file(std::size_t core_index = 0) const noexcept;
  [[nodiscard]] AddressSpace& address_space() noexcept;
  [[nodiscard]] const AddressSpace& address_space() const noexcept;
  [[nodiscard]] SyscallResult dispatch_syscall(Context&, SyscallKind) noexcept;
  [[nodiscard]] CpuidResult dispatch_cpuid(Context&, CpuidRequest) noexcept;
  void set_pending_stop(RunResult);

private:
  HostCallbacks& callbacks_;
  Options options_;
  std::unique_ptr<AddressSpace> address_space_;
  std::unique_ptr<MachineImpl> machine_;
  std::optional<RunResult> pending_stop_;
};

} // namespace x86sim

template<>
struct std::formatter<x86sim::Protection> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::Protection protection, FormatContext& ctx) const {
    const auto bits = static_cast<unsigned>(protection);
    if ((bits & ~0x7u) != 0)
      return std::format_to(ctx.out(), "unknown(0x{:02x})", bits);
    if (protection == x86sim::Protection::none)
      return std::format_to(ctx.out(), "none");

    auto out = ctx.out();
    if (x86sim::has_protection(protection, x86sim::Protection::read))
      out = std::format_to(out, "r");
    if (x86sim::has_protection(protection, x86sim::Protection::write))
      out = std::format_to(out, "w");
    if (x86sim::has_protection(protection, x86sim::Protection::execute))
      out = std::format_to(out, "x");
    return out;
  }
};

template<>
struct std::formatter<x86sim::CoreModel> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::CoreModel model, FormatContext& ctx) const {
    using enum x86sim::CoreModel;
    switch (model) {
    case out_of_order:
      return std::format_to(ctx.out(), "out_of_order");
    case sequential:
      return std::format_to(ctx.out(), "sequential");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(model));
  }
};

template<>
struct std::formatter<x86sim::SyscallKind> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::SyscallKind kind, FormatContext& ctx) const {
    using enum x86sim::SyscallKind;
    switch (kind) {
    case int80:
      return std::format_to(ctx.out(), "int80");
    case syscall64:
      return std::format_to(ctx.out(), "syscall64");
    case sysenter:
      return std::format_to(ctx.out(), "sysenter");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(kind));
  }
};

template<>
struct std::formatter<x86sim::StopReason> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::StopReason reason, FormatContext& ctx) const {
    using enum x86sim::StopReason;
    switch (reason) {
    case guest_exit:
      return std::format_to(ctx.out(), "guest_exit");
    case instruction_limit:
      return std::format_to(ctx.out(), "instruction_limit");
    case x86_exception:
      return std::format_to(ctx.out(), "x86_exception");
    case host_request:
      return std::format_to(ctx.out(), "host_request");
    case unsupported_syscall:
      return std::format_to(ctx.out(), "unsupported_syscall");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(reason));
  }
};

template<>
struct std::formatter<x86sim::MemoryError> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::MemoryError error, FormatContext& ctx) const {
    using enum x86sim::MemoryError;
    switch (error) {
    case unaligned_address:
      return std::format_to(ctx.out(), "address is not page-aligned");
    case zero_size:
      return std::format_to(ctx.out(), "range has zero size");
    case unmapped_address:
      return std::format_to(ctx.out(), "address is not mapped");
    case out_of_memory:
      return std::format_to(ctx.out(), "out of memory");
    case mapping_failed:
      return std::format_to(ctx.out(), "mapping failed");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(error));
  }
};

template<>
struct std::formatter<x86sim::CpuidRequest> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::CpuidRequest& request, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "function=0x{:08x}, subfunction=0x{:08x}", request.function, request.subfunction);
  }
};

template<>
struct std::formatter<x86sim::CpuidResult> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::CpuidResult& result, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "eax=0x{:08x}, ebx=0x{:08x}, ecx=0x{:08x}, edx=0x{:08x}", result.eax,
                          result.ebx, result.ecx, result.edx);
  }
};

template<>
struct std::formatter<x86sim::SyscallResult> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::SyscallResult& result, FormatContext& ctx) const {
    auto out = std::format_to(ctx.out(), "reason={}, continue_execution={}", result.reason,
                              result.continue_execution);
    if (!result.message.empty())
      out = std::format_to(out, ", message={}", result.message);
    return out;
  }
};

template<>
struct std::formatter<x86sim::Options> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::Options& options, FormatContext& ctx) const {
    return std::format_to(ctx.out(),
                          "core={}, core_count={}, sse={}, x87={}, perfect_cache={}, static_branchpred={}, "
                          "log_filename={}",
                          options.core, options.core_count, options.sse, options.x87, options.perfect_cache,
                          options.static_branchpred, options.log_filename.string());
  }
};

template<>
struct std::formatter<x86sim::RunOptions> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::RunOptions& options, FormatContext& ctx) const {
    if (options.instruction_limit)
      return std::format_to(ctx.out(), "instruction_limit={}", *options.instruction_limit);
    return std::format_to(ctx.out(), "instruction_limit=none");
  }
};

template<>
struct std::formatter<x86sim::Stats> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::Stats& stats, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "{} cycles, {} instructions", stats.cycles, stats.instructions);
  }
};

template<>
struct std::formatter<x86sim::X86Exception> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::X86Exception& exception, FormatContext& ctx) const {
    auto out = ctx.out();
    if (!exception.message.empty()) {
      out = std::format_to(out, "{}", exception.message);
    } else {
      out = std::format_to(out, "Exception {} code={} addr=0x{:x} @ rip 0x{:x}",
                           static_cast<unsigned>(exception.vector), exception.error_code,
                           exception.virtual_address, exception.rip);
    }
    if (!exception.context.empty())
      out = std::format_to(out, "\n{}", exception.context);
    return out;
  }
};

template<>
struct std::formatter<x86sim::RunResult> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::RunResult& result, FormatContext& ctx) const {
    auto out = std::format_to(ctx.out(), "{} after {}", result.reason, result.stats);
    if (result.x86_exception)
      return std::format_to(out, "\n{}", *result.x86_exception);
    if (!result.message.empty())
      out = std::format_to(out, ": {}", result.message);
    return out;
  }
};

template<>
struct std::formatter<x86sim::Machine> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::Machine& machine, FormatContext& ctx) const {
    const auto& options = machine.options();
    const auto stats = machine.stats();
    return std::format_to(ctx.out(),
                          "x86sim::Machine(core={}, core_count={}, address_space={}, cycles={}, instructions={}, "
                          "sse={}, x87={}, perfect_cache={}, static_branch_prediction={}, core0={})",
                          options.core, options.core_count, static_cast<const void*>(&machine.address_space()), stats.cycles,
                          stats.instructions, options.sse, options.x87, options.perfect_cache,
                          options.static_branchpred,
                          options.core_count == 0 ? x86sim::RegisterFile{} : machine.register_file(0));
  }
};

#endif
