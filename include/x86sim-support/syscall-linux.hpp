#ifndef X86SIM_SUPPORT_SYSCALL_LINUX_HPP
#define X86SIM_SUPPORT_SYSCALL_LINUX_HPP

#include "x86sim-support/cpuid.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace x86sim::linux_syscalls {

// Identity the support library assigns to each guest process/thread it manages.
// This is owned by the support library, not the Machine: the caller keeps a
// CpuState + AddressSpace per process and tags them with a ProcessId so the
// process-aware syscalls (fork/wait/exit/getpid, per-process fd tables, ...)
// can find their bookkeeping.
using ProcessId = std::uint64_t;
inline constexpr ProcessId invalid_process_id = std::numeric_limits<ProcessId>::max();
inline constexpr ProcessId initial_process_id = 1;

struct ProcessTable;

[[nodiscard]] std::shared_ptr<ProcessTable> make_process_table();

namespace detail {

inline constexpr address_t default_mmap_base = 0x700000000000ULL;
inline constexpr address_t default_brk_base = 0x600000000000ULL;

[[nodiscard]] SyscallResult unsupported_syscall(const CpuState&, SyscallKind);

} // namespace detail

// Public, minimal ABI helpers for callers that drive the syscall instruction
// themselves (e.g. routing read/write/exit to host-supplied streams) without
// pulling in the internal detail surface. Thin forwarders over detail::.
namespace abi {

[[nodiscard]] word_t syscall_number(const RegisterFile& context) noexcept;
[[nodiscard]] word_t syscall_arg(const CpuState& context, std::size_t index) noexcept;
[[nodiscard]] SyscallResult return_value(CpuState& context, std::int64_t value) noexcept;

inline constexpr word_t read = 0;
inline constexpr word_t write = 1;
inline constexpr word_t exit = 60;
inline constexpr word_t exit_group = 231;

} // namespace abi

template<typename Handler>
concept SyscallHandler = requires(Handler& handler, Machine& machine, ProcessId pid, CpuState& context,
                                  AddressSpace& space, SyscallKind kind) {
  { handler.try_syscall(machine, pid, context, space, kind) } -> std::same_as<std::optional<SyscallResult>>;
};

struct SysProcess {
  explicit SysProcess(std::shared_ptr<ProcessTable> processes = make_process_table()) noexcept
      : processes(std::move(processes)) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  std::shared_ptr<ProcessTable> processes;
};

struct SysFutex {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysSignals {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysBrk {
  constexpr explicit SysBrk(address_t initial_break = detail::default_brk_base) noexcept
      : initial_break(initial_break) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  address_t initial_break;
};

struct SysArchPrctl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysGetIdentity {
  explicit SysGetIdentity(std::shared_ptr<ProcessTable> processes = make_process_table()) noexcept
      : processes(std::move(processes)) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  std::shared_ptr<ProcessTable> processes;
};

struct SysUname {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysGetcwd {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysSetTidAddress {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysSetRobustList {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysRseq {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

// getrandom(buf, buflen, flags): fills the guest buffer with deterministic
// pseudo-random bytes (a fixed-seed splitmix64 stream) and returns buflen. The
// fill is portable and reproducible -- no host entropy source is touched -- so
// the simulator stays deterministic. glibc >= 2.36 calls this at startup to seed
// the stack canary / pointer guard.
struct SysGetrandom {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysPrlimit64 {
  struct Limit {
    word_t current;
    word_t maximum;
  };

  constexpr SysPrlimit64() noexcept
      : limits{Limit{unlimited, unlimited},
               Limit{unlimited, unlimited},
               Limit{unlimited, unlimited},
               Limit{8 * 1024 * 1024, unlimited},
               Limit{0, unlimited},
               Limit{unlimited, unlimited},
               Limit{4096, 4096},
               Limit{1024, 1024},
               Limit{64 * 1024, 64 * 1024},
               Limit{unlimited, unlimited},
               Limit{unlimited, unlimited},
               Limit{4096, 4096},
               Limit{819200, 819200},
               Limit{0, 0},
               Limit{0, 0},
               Limit{unlimited, unlimited}} {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  static constexpr word_t unlimited = ~word_t{0};
  std::array<Limit, 16> limits;
};

struct SysExit {
  explicit SysExit(int* exit_status = nullptr, std::shared_ptr<ProcessTable> processes = make_process_table()) noexcept
      : exit_status(exit_status), processes(std::move(processes)) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  int* exit_status;
  std::shared_ptr<ProcessTable> processes;
};

struct SysExitGroup {
  explicit SysExitGroup(int* exit_status = nullptr,
                        std::shared_ptr<ProcessTable> processes = make_process_table()) noexcept
      : exit_status(exit_status), processes(std::move(processes)) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  int* exit_status;
  std::shared_ptr<ProcessTable> processes;
};

struct SysMmap {
  constexpr explicit SysMmap(address_t next_mapping_address = detail::default_mmap_base) noexcept
      : next_mapping_address(next_mapping_address) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  address_t next_mapping_address;
};

struct SysMunmap {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysMprotect {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysMremap {
  constexpr explicit SysMremap(address_t next_mapping_address = detail::default_mmap_base + 0x8000000000ULL) noexcept
      : next_mapping_address(next_mapping_address) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  address_t next_mapping_address;
};

template<typename... Handlers>
class Chain {
public:
  constexpr Chain() = default;

  constexpr explicit Chain(Handlers... handlers) : handlers_(std::move(handlers)...) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine& machine, ProcessId pid, CpuState& context,
                                                         AddressSpace& space, SyscallKind kind) noexcept {
    return try_syscall_at<0>(machine, pid, context, space, kind);
  }

  [[nodiscard]] SyscallResult operator()(Machine& machine, ProcessId pid, CpuState& context, AddressSpace& space,
                                         SyscallKind kind) {
    if (auto result = try_syscall(machine, pid, context, space, kind))
      return *result;
    return detail::unsupported_syscall(context, kind);
  }

  [[nodiscard]] constexpr auto into_tuple() const& { return handlers_; }

  [[nodiscard]] constexpr auto into_tuple() && { return std::move(handlers_); }

private:
  template<std::size_t index>
  [[nodiscard]] std::optional<SyscallResult> try_syscall_at(Machine& machine, ProcessId pid, CpuState& context,
                                                            AddressSpace& space, SyscallKind kind) noexcept {
    if constexpr (index == sizeof...(Handlers)) {
      return std::nullopt;
    } else {
      if (auto result = std::get<index>(handlers_).try_syscall(machine, pid, context, space, kind))
        return result;
      return try_syscall_at<index + 1>(machine, pid, context, space, kind);
    }
  }

  std::tuple<Handlers...> handlers_;
};

namespace detail {

template<typename T>
[[nodiscard]] constexpr auto into_tuple(T&& value) {
  if constexpr (requires { std::forward<T>(value).into_tuple(); })
    return std::forward<T>(value).into_tuple();
  else
    return std::tuple<std::decay_t<T>>(std::forward<T>(value));
}

template<typename Tuple, std::size_t... indexes>
[[nodiscard]] constexpr auto chain_from_tuple(Tuple&& tuple, std::index_sequence<indexes...>) {
  using TupleType = std::remove_cvref_t<Tuple>;
  return Chain<std::tuple_element_t<indexes, TupleType>...>(std::get<indexes>(std::forward<Tuple>(tuple))...);
}

template<typename Tuple>
[[nodiscard]] constexpr auto chain_from_tuple(Tuple&& tuple) {
  using TupleType = std::remove_cvref_t<Tuple>;
  return chain_from_tuple(std::forward<Tuple>(tuple), std::make_index_sequence<std::tuple_size_v<TupleType>>{});
}

} // namespace detail

template<SyscallHandler Lhs, SyscallHandler Rhs>
[[nodiscard]] constexpr auto operator|(Lhs&& lhs, Rhs&& rhs) {
  return detail::chain_from_tuple(
      std::tuple_cat(detail::into_tuple(std::forward<Lhs>(lhs)), detail::into_tuple(std::forward<Rhs>(rhs))));
}

template<SyscallHandler Syscalls>
class Host : public HostCallbacks {
public:
  constexpr explicit Host(Syscalls syscalls = {}) : syscalls_(std::move(syscalls)) {}

  [[nodiscard]] CpuidResult cpuid(Machine&, CpuState&, AddressSpace&, CpuidRequest request) noexcept override {
    return defaults::default_cpuid(request);
  }

  [[nodiscard]] SyscallResult syscall(Machine& machine, CpuState& context, AddressSpace& space,
                                      SyscallKind kind) override {
    // Phase 1 runs a single process; Phase 2's scheduler will supply the live
    // ProcessId here instead of the fixed initial id.
    return syscalls_(machine, initial_process_id, context, space, kind);
  }

private:
  Syscalls syscalls_;
};

template<SyscallHandler Syscalls>
Host(Syscalls) -> Host<std::decay_t<Syscalls>>;

template<SyscallHandler Syscalls>
[[nodiscard]] constexpr auto host(Syscalls&& syscalls) {
  return Host<std::decay_t<Syscalls>>(std::forward<Syscalls>(syscalls));
}

} // namespace x86sim::linux_syscalls

#endif
