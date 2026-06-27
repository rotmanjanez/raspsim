#ifndef X86SIM_SUPPORT_SYSCALL_LINUX_POSIX_HPP
#define X86SIM_SUPPORT_SYSCALL_LINUX_POSIX_HPP

// Host-POSIX Linux syscall handlers: file, directory, socket, and host-time I/O
// backed by the host operating system. These live in the x86sim::posix library
// (separate from the portable x86sim::defaults library). Callers that only need
// the portable handlers (decode helpers, heap emulation, register-only syscalls)
// should depend on syscall-linux.hpp / x86sim::defaults instead.
#include "x86sim-support/syscall-linux.hpp"

#include <optional>

namespace x86sim::linux_syscalls {

struct SysRead {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysReadlink {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysWrite {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysPreadPwrite {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysOpen {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysClose {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysPipe {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysFstat {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysStat {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysFileSystem {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysLseek {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysIoctl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysFcntl {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysSocket {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysConnect {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysSelect {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysGetdents64 {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

struct SysTime {
  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;
};

} // namespace x86sim::linux_syscalls

#endif
