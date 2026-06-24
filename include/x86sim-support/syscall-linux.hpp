#ifndef X86SIM_SUPPORT_SYSCALL_LINUX_HPP
#define X86SIM_SUPPORT_SYSCALL_LINUX_HPP

#include "x86sim-support/defaults.hpp"

#include <array>
#include <concepts>
#include <cstdint>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace x86sim::linux_syscalls {

namespace detail {

inline constexpr address_t default_mmap_base = 0x700000000000ULL;
inline constexpr address_t default_brk_base = 0x600000000000ULL;

[[nodiscard]] SyscallResult unsupported_syscall(const RegisterFile&, SyscallKind);

} // namespace detail

template<typename Handler>
concept SyscallHandler = requires(Handler& handler, Machine& machine, RegisterFile& context, SyscallKind kind) {
  { handler.try_syscall(machine, context, kind) } -> std::same_as<std::optional<SyscallResult>>;
};

struct SysRead {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysReadlink {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysWrite {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysPreadPwrite {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysOpen {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysClose {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysFstat {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysStat {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysFileSystem {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysLseek {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysIoctl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysFcntl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysSocket {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysConnect {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysSignals {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysBrk {
  constexpr explicit SysBrk(address_t initial_break = detail::default_brk_base) noexcept
      : minimum_break(initial_break), current_break(initial_break), mapped_end(initial_break) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;

  address_t minimum_break;
  address_t current_break;
  address_t mapped_end;
};

struct SysArchPrctl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysGetIdentity {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysUname {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysGetcwd {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysSetTidAddress {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysSetRobustList {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysRseq {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

struct SysPrlimit64 {
  struct Limit {
    word_t current;
    word_t maximum;
  };

  constexpr SysPrlimit64() noexcept
      : limits{Limit{unlimited, unlimited}, Limit{unlimited, unlimited}, Limit{unlimited, unlimited},
               Limit{8 * 1024 * 1024, unlimited}, Limit{0, unlimited}, Limit{unlimited, unlimited},
               Limit{4096, 4096}, Limit{1024, 1024}, Limit{64 * 1024, 64 * 1024}, Limit{unlimited, unlimited},
               Limit{unlimited, unlimited}, Limit{4096, 4096}, Limit{819200, 819200}, Limit{0, 0}, Limit{0, 0},
               Limit{unlimited, unlimited}} {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;

  static constexpr word_t unlimited = ~word_t{0};
  std::array<Limit, 16> limits;
};

struct SysExit {
  constexpr explicit SysExit(int* exit_status = nullptr) noexcept : exit_status(exit_status) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;

  int* exit_status;
};

struct SysExitGroup {
  constexpr explicit SysExitGroup(int* exit_status = nullptr) noexcept : exit_status(exit_status) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;

  int* exit_status;
};

struct SysMmap {
  constexpr explicit SysMmap(address_t next_mapping_address = detail::default_mmap_base) noexcept
      : next_mapping_address(next_mapping_address) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;

  address_t next_mapping_address;
};

struct SysMunmap {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, RegisterFile&, SyscallKind) noexcept;
};

template<typename... Handlers>
class Chain {
public:
  constexpr Chain() = default;

  constexpr explicit Chain(Handlers... handlers) : handlers_(std::move(handlers)...) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
    return try_syscall_at<0>(machine, context, kind);
  }

  [[nodiscard]] SyscallResult operator()(Machine& machine, RegisterFile& context, SyscallKind kind) {
    if (auto result = try_syscall(machine, context, kind))
      return *result;
    return detail::unsupported_syscall(context, kind);
  }

  [[nodiscard]] constexpr auto into_tuple() const& {
    return handlers_;
  }

  [[nodiscard]] constexpr auto into_tuple() && {
    return std::move(handlers_);
  }

private:
  template<std::size_t index>
  [[nodiscard]] std::optional<SyscallResult> try_syscall_at(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
    if constexpr (index == sizeof...(Handlers)) {
      return std::nullopt;
    } else {
      if (auto result = std::get<index>(handlers_).try_syscall(machine, context, kind))
        return result;
      return try_syscall_at<index + 1>(machine, context, kind);
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

  [[nodiscard]] CpuidResult cpuid(Machine&, RegisterFile&, CpuidRequest request) noexcept override {
    return defaults::default_cpuid(request);
  }

  [[nodiscard]] SyscallResult syscall(Machine& machine, RegisterFile& context, SyscallKind kind) override {
    return syscalls_(machine, context, kind);
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
