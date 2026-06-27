#ifndef RASPSIM_RASPSIM_HPP
#define RASPSIM_RASPSIM_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <string>

#include "x86sim/addrspace.hpp"
#include "x86sim/registerfile.hpp"

namespace x86sim {

struct MachineImpl;
struct CpuState;

enum class CoreModel { out_of_order, sequential };
enum class SegmentRegister : std::uint8_t { es = 0, cs = 1, ss = 2, ds = 3, fs = 4, gs = 5 };

inline constexpr std::size_t segment_register_count = 6;

[[nodiscard]] constexpr std::size_t segment_register_index(SegmentRegister reg) noexcept {
  return static_cast<std::size_t>(reg);
}

enum class SyscallKind { int80, syscall64, sysenter };
enum class StopReason { guest_exit, instruction_limit, x86_exception, host_request, unsupported_syscall };

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
  virtual SyscallResult syscall(Machine&, CpuState&, AddressSpace&, SyscallKind) = 0;
  virtual CpuidResult cpuid(Machine&, CpuState&, AddressSpace&, CpuidRequest) noexcept = 0;
};

struct LogOptions {
  static constexpr address_t invalid_rip = std::numeric_limits<address_t>::max();

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
};

struct DebugOptions {
  static constexpr address_t invalid_rip = std::numeric_limits<address_t>::max();

  bool dump_state_now = false;
  bool abort_at_end = false;

  address_t start_at_rip = 0;
  bool include_dyn_linker = false;
  bool trigger_mode = false;
  std::uint64_t pause_at_startup = 0;

  bool perfect_cache = false;
  bool static_branchpred = false;

  std::string bbcache_dump_filename;
  std::uint64_t sequential_mode_insns = 0;
  bool exit_after_fullsim = false;
};

struct Options {
  static constexpr address_t invalid_rip = LogOptions::invalid_rip;

  CoreModel core = CoreModel::out_of_order;
  std::uint64_t core_count = 1;
  bool sse = true;
  bool x87 = true;
  bool quiet = false;

  std::uint64_t stop_at_user_insns = std::numeric_limits<std::uint64_t>::max();
  std::uint64_t stop_at_iteration = std::numeric_limits<std::uint64_t>::max();
  address_t stop_at_rip = invalid_rip;

  LogOptions log;
  DebugOptions debug;
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

// The complete externally-visible architectural state of a single CPU, and the
// single definition of that state: the simulator's internal context derives
// from it, so loading/storing it across a run is a plain slice/assign rather
// than a field-by-field copy. Inherits RegisterFile, so the general-purpose/SSE
// register accessors (operator[]) are available directly. The caller owns one
// or more of these (one per process/thread) and hands them to Machine::run.
struct CpuState : RegisterFile {
  std::array<address_t, segment_register_count> segment_bases{};
  std::array<std::uint16_t, segment_register_count> segment_selectors{};
  address_t swapgs_base = 0;
  std::array<word_t, 8> fpstack{};
  word_t fpcw = 0;
  std::uint32_t mxcsr = 0;
  bool use32 = true;
  bool use64 = true;
  address_t virt_addr_mask = std::numeric_limits<address_t>::max();
  std::uint32_t internal_eflags = 0;
};

class Machine {
public:
  explicit Machine(HostCallbacks&, Options = {});
  ~Machine();

  // Execute `state` against `space` until a stop is reached, mutating `state`
  // in place and returning why execution stopped. The caller owns the CpuState
  // and AddressSpace and may keep and swap several of them to model multiple
  // processes/threads on the same Machine.
  [[nodiscard]] RunResult run(CpuState& state, AddressSpace& space, RunOptions = {});

  [[nodiscard]] const Options& options() const noexcept;
  [[nodiscard]] Stats stats() const noexcept;

private:
  // The core delivers guest syscalls/cpuid to callbacks_ and records the
  // resulting stop in MachineImpl::pending_stop, so it owns those internals.
  friend struct MachineImpl;
  // An AddressSpace bound to this Machine reads machine_ to reach the code cache
  // it must invalidate on map/unmap/protect.
  friend class AddressSpace;

  HostCallbacks& callbacks_;
  Options options_;
  std::unique_ptr<MachineImpl> machine_;
};

} // namespace x86sim

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
struct std::formatter<x86sim::SegmentRegister> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::SegmentRegister reg, FormatContext& ctx) const {
    using enum x86sim::SegmentRegister;
    switch (reg) {
    case es:
      return std::format_to(ctx.out(), "es");
    case cs:
      return std::format_to(ctx.out(), "cs");
    case ss:
      return std::format_to(ctx.out(), "ss");
    case ds:
      return std::format_to(ctx.out(), "ds");
    case fs:
      return std::format_to(ctx.out(), "fs");
    case gs:
      return std::format_to(ctx.out(), "gs");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(reg));
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
    return std::format_to(ctx.out(), "eax=0x{:08x}, ebx=0x{:08x}, ecx=0x{:08x}, edx=0x{:08x}", result.eax, result.ebx,
                          result.ecx, result.edx);
  }
};

template<>
struct std::formatter<x86sim::SyscallResult> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::SyscallResult& result, FormatContext& ctx) const {
    auto out = std::format_to(ctx.out(), "reason={}, continue_execution={}", result.reason, result.continue_execution);
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
                          options.core, options.core_count, options.sse, options.x87, options.debug.perfect_cache,
                          options.debug.static_branchpred, options.log.log_filename.string());
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
      out =
          std::format_to(out, "Exception {} code={} addr=0x{:x} @ rip 0x{:x}", static_cast<unsigned>(exception.vector),
                         exception.error_code, exception.virtual_address, exception.rip);
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
                          "x86sim::Machine(core={}, core_count={}, cycles={}, instructions={}, sse={}, x87={}, "
                          "perfect_cache={}, static_branch_prediction={})",
                          options.core, options.core_count, stats.cycles, stats.instructions, options.sse, options.x87,
                          options.debug.perfect_cache, options.debug.static_branchpred);
  }
};

#endif
