// Internal shared surface for the Linux syscall support library.
//
// This header is NOT installed. It is included by both translation units that
// make up the support library:
//   * syscall-linux.cpp       — the fully portable handlers + helper definitions
//   * syscall-linux-posix.cpp — the host-POSIX handlers (file/socket/time I/O)
//
// It holds the shared decode/error helpers, the Linux ABI constants, the
// ProcessTable definition, and the brk bookkeeping. Functions declared here are
// in namespace detail with external linkage; their single definitions live in
// the portable translation unit so the POSIX unit can reuse them by linking the
// portable library. Anything that needs a host header (e.g. host_errno_to_linux)
// stays in the POSIX unit and is deliberately absent here.
#ifndef X86SIM_SUPPORT_SYSCALL_LINUX_DETAIL_HPP
#define X86SIM_SUPPORT_SYSCALL_LINUX_DETAIL_HPP

#include "x86sim-support/syscall-linux.hpp"
#include "x86sim/addrspace.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace x86sim::linux_syscalls {

struct ProcessTable {
  struct ExitRecord {
    ProcessId child = invalid_process_id;
    int status = 0;
  };

  std::vector<std::unique_ptr<AddressSpace>> address_spaces;
  std::unordered_map<ProcessId, ProcessId> parent;
  std::unordered_map<ProcessId, ProcessId> vfork_parent;
  std::unordered_map<ProcessId, std::vector<ProcessId>> live_children;
  std::unordered_map<ProcessId, std::vector<ExitRecord>> exited_children;

  struct WaitRequest {
    std::int64_t requested_pid = -1;
    address_t status_address = 0;
  };
  std::unordered_map<ProcessId, WaitRequest> waiters;
};

namespace detail {

inline constexpr word_t syscall_read = 0;
inline constexpr word_t syscall_write = 1;
inline constexpr word_t syscall_open = 2;
inline constexpr word_t syscall_close = 3;
inline constexpr word_t syscall_stat = 4;
inline constexpr word_t syscall_fstat = 5;
inline constexpr word_t syscall_lstat = 6;
inline constexpr word_t syscall_lseek = 8;
inline constexpr word_t syscall_mmap = 9;
inline constexpr word_t syscall_munmap = 11;
inline constexpr word_t syscall_brk = 12;
inline constexpr word_t syscall_rt_sigaction = 13;
inline constexpr word_t syscall_rt_sigprocmask = 14;
inline constexpr word_t syscall_ioctl = 16;
inline constexpr word_t syscall_pread64 = 17;
inline constexpr word_t syscall_pwrite64 = 18;
inline constexpr word_t syscall_access = 21;
inline constexpr word_t syscall_pipe = 22;
inline constexpr word_t syscall_select = 23;
inline constexpr word_t syscall_mremap = 25;
inline constexpr word_t syscall_dup = 32;
inline constexpr word_t syscall_dup2 = 33;
inline constexpr word_t syscall_nanosleep = 35;
inline constexpr word_t syscall_clone = 56;
inline constexpr word_t syscall_fork = 57;
inline constexpr word_t syscall_vfork = 58;
inline constexpr word_t syscall_execve = 59;
inline constexpr word_t syscall_wait4 = 61;
inline constexpr word_t syscall_fcntl = 72;
inline constexpr word_t syscall_truncate = 76;
inline constexpr word_t syscall_ftruncate = 77;
inline constexpr word_t syscall_chdir = 80;
inline constexpr word_t syscall_fchdir = 81;
inline constexpr word_t syscall_rename = 82;
inline constexpr word_t syscall_mkdir = 83;
inline constexpr word_t syscall_rmdir = 84;
inline constexpr word_t syscall_unlink = 87;
inline constexpr word_t syscall_umask = 95;
inline constexpr word_t syscall_gettimeofday = 96;
inline constexpr word_t syscall_getpid = 39;
inline constexpr word_t syscall_socket = 41;
inline constexpr word_t syscall_connect = 42;
inline constexpr word_t syscall_exit = 60;
inline constexpr word_t syscall_uname = 63;
inline constexpr word_t syscall_getcwd = 79;
inline constexpr word_t syscall_readlink = 89;
inline constexpr word_t syscall_chmod = 90;
inline constexpr word_t syscall_fchmod = 91;
inline constexpr word_t syscall_getrusage = 98;
inline constexpr word_t syscall_getuid = 102;
inline constexpr word_t syscall_getgid = 104;
inline constexpr word_t syscall_geteuid = 107;
inline constexpr word_t syscall_getegid = 108;
inline constexpr word_t syscall_utime = 132;
inline constexpr word_t syscall_arch_prctl = 158;
inline constexpr word_t syscall_gettid = 186;
inline constexpr word_t syscall_futex = 202;
inline constexpr word_t syscall_getdents64 = 217;
inline constexpr word_t syscall_set_tid_address = 218;
inline constexpr word_t syscall_clock_gettime = 228;
inline constexpr word_t syscall_exit_group = 231;
inline constexpr word_t syscall_openat = 257;
inline constexpr word_t syscall_mkdirat = 258;
inline constexpr word_t syscall_newfstatat = 262;
inline constexpr word_t syscall_unlinkat = 263;
inline constexpr word_t syscall_renameat = 264;
inline constexpr word_t syscall_readlinkat = 267;
inline constexpr word_t syscall_faccessat = 269;
inline constexpr word_t syscall_set_robust_list = 273;
inline constexpr word_t syscall_prlimit64 = 302;
inline constexpr word_t syscall_dup3 = 292;
inline constexpr word_t syscall_rseq = 334;

inline constexpr int linux_enoent = 2;
inline constexpr int linux_esrch = 3;
inline constexpr int linux_eio = 5;
inline constexpr int linux_e2big = 7;
inline constexpr int linux_enoexec = 8;
inline constexpr int linux_ebadf = 9;
inline constexpr int linux_echild = 10;
inline constexpr int linux_eagain = 11;
inline constexpr int linux_enomem = 12;
inline constexpr int linux_eacces = 13;
inline constexpr int linux_efault = 14;
inline constexpr int linux_enotty = 25;
inline constexpr int linux_einval = 22;
inline constexpr int linux_emfile = 24;
inline constexpr int linux_eafnosupport = 97;
inline constexpr int linux_enametoolong = 36;
inline constexpr int linux_enosys = 38;
inline constexpr int linux_etimedout = 110;

inline constexpr word_t linux_o_accmode = 00000003;
inline constexpr word_t linux_o_rdonly = 00000000;
inline constexpr word_t linux_o_wronly = 00000001;
inline constexpr word_t linux_o_rdwr = 00000002;
inline constexpr word_t linux_o_creat = 00000100;
inline constexpr word_t linux_o_excl = 00000200;
inline constexpr word_t linux_o_noctty = 00000400;
inline constexpr word_t linux_o_trunc = 00001000;
inline constexpr word_t linux_o_append = 00002000;
inline constexpr word_t linux_o_nonblock = 00004000;
inline constexpr word_t linux_o_dsync = 00010000;
inline constexpr word_t linux_o_largefile = 00100000;
inline constexpr word_t linux_o_directory = 00200000;
inline constexpr word_t linux_o_nofollow = 00400000;
inline constexpr word_t linux_o_cloexec = 02000000;
inline constexpr word_t linux_o_sync = 04010000;

inline constexpr word_t linux_at_fdcwd = static_cast<word_t>(-100);
inline constexpr word_t linux_at_symlink_nofollow = 0x100;
inline constexpr word_t linux_at_removedir = 0x200;
inline constexpr word_t linux_at_empty_path = 0x1000;

inline constexpr word_t linux_f_dupfd = 0;
inline constexpr word_t linux_f_getfd = 1;
inline constexpr word_t linux_f_setfd = 2;
inline constexpr word_t linux_f_getfl = 3;
inline constexpr word_t linux_f_setfl = 4;
inline constexpr word_t linux_f_getlk = 5;
inline constexpr word_t linux_f_setlk = 6;
inline constexpr word_t linux_f_setlkw = 7;
inline constexpr word_t linux_f_dupfd_cloexec = 1030;

inline constexpr word_t linux_fd_cloexec = 1;
inline constexpr word_t linux_f_rdlck = 0;
inline constexpr word_t linux_f_wrlck = 1;
inline constexpr word_t linux_f_unlck = 2;

inline constexpr word_t linux_prot_read = 0x1;
inline constexpr word_t linux_prot_write = 0x2;
inline constexpr word_t linux_prot_exec = 0x4;

inline constexpr word_t linux_map_shared = 0x01;
inline constexpr word_t linux_map_private = 0x02;
inline constexpr word_t linux_map_fixed = 0x10;
inline constexpr word_t linux_map_anonymous = 0x20;
inline constexpr word_t linux_map_denied_write = 0x800;
inline constexpr word_t linux_map_executable = 0x1000;
inline constexpr word_t linux_map_noreserve = 0x4000;
inline constexpr word_t linux_map_populate = 0x8000;
inline constexpr word_t linux_map_stack = 0x20000;

inline constexpr word_t linux_mremap_maymove = 1;
inline constexpr word_t linux_mremap_fixed = 2;
inline constexpr word_t linux_mremap_dontunmap = 4;

inline constexpr word_t linux_af_unix = 1;
inline constexpr word_t linux_af_inet = 2;
inline constexpr word_t linux_af_inet6 = 10;

inline constexpr word_t linux_sock_type_mask = 0xf;
inline constexpr word_t linux_sock_stream = 1;
inline constexpr word_t linux_sock_dgram = 2;
inline constexpr word_t linux_sock_raw = 3;
inline constexpr word_t linux_sock_seqpacket = 5;
inline constexpr word_t linux_sock_nonblock = 00004000;
inline constexpr word_t linux_sock_cloexec = 02000000;

inline constexpr word_t linux_arch_set_gs = 0x1001;
inline constexpr word_t linux_arch_set_fs = 0x1002;
inline constexpr word_t linux_arch_get_fs = 0x1003;
inline constexpr word_t linux_arch_get_gs = 0x1004;

inline constexpr word_t linux_futex_wait = 0;
inline constexpr word_t linux_futex_wake = 1;
inline constexpr word_t linux_futex_wait_bitset = 9;
inline constexpr word_t linux_futex_wake_bitset = 10;
inline constexpr word_t linux_futex_private_flag = 128;
inline constexpr word_t linux_futex_clock_realtime = 256;
inline constexpr word_t linux_futex_cmd_mask = ~(linux_futex_private_flag | linux_futex_clock_realtime);

inline constexpr word_t synthetic_pid = 1;
inline constexpr word_t synthetic_uid = 1000;
inline constexpr word_t synthetic_gid = 1000;

inline constexpr std::size_t max_path_length = 4096;
inline constexpr std::size_t io_chunk_size = 64 * 1024;

// --- shared decode / error helpers (defined in the portable TU) -----------
[[nodiscard]] word_t syscall_number(const RegisterFile& context) noexcept;
[[nodiscard]] word_t syscall_arg(const CpuState& context, std::size_t index) noexcept;
[[nodiscard]] bool handles(const CpuState& context, SyscallKind kind, word_t number) noexcept;
[[nodiscard]] SyscallResult return_value(CpuState& context, std::int64_t value);
[[nodiscard]] SyscallResult return_error(CpuState& context, int error);
[[nodiscard]] int memory_error_to_linux(MemoryError error) noexcept;
[[nodiscard]] std::optional<int> checked_fd(word_t raw_fd) noexcept;
[[nodiscard]] bool fits_host_transfer(word_t count) noexcept;
[[nodiscard]] bool range_overflows(address_t start, word_t length) noexcept;
[[nodiscard]] bool is_linux_at_fdcwd(word_t raw_fd) noexcept;
[[nodiscard]] std::uint64_t read_le(std::span<const std::byte> bytes, std::size_t offset, std::size_t width) noexcept;

struct CStringResult {
  bool ok = false;
  int error = 0;
  std::string value;
};

[[nodiscard]] CStringResult read_c_string(AddressSpace& space, address_t address,
                                          std::size_t max_length = max_path_length) noexcept;

template<std::size_t Size>
[[nodiscard]] std::optional<int> read_guest_memory(AddressSpace& space, address_t address,
                                                   std::array<std::byte, Size>& out) noexcept {
  auto read = space.read(address, out);
  if (!read)
    return memory_error_to_linux(read.error());
  return std::nullopt;
}

template<std::size_t Size>
void write_le(std::array<std::byte, Size>& bytes, std::size_t offset, std::uint64_t value, std::size_t width) noexcept {
  for (std::size_t i = 0; i < width; ++i)
    bytes[offset + i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);
}

// --- brk bookkeeping ------------------------------------------------------
// Per-process heap break state shared by SysBrk and the process lifecycle
// helpers. Declarations live here; the map and definitions live in the portable
// TU.
struct BrkState {
  address_t minimum_break = default_brk_base;
  address_t current_break = default_brk_base;
  address_t mapped_end = default_brk_base;
};

[[nodiscard]] BrkState& brk_state(ProcessId context_id, address_t initial_break = default_brk_base);
void clone_brk_state(ProcessId parent, ProcessId child);
void reset_brk_state(ProcessId context_id, address_t initial_break = default_brk_base);
void erase_brk_state(ProcessId context_id);

// --- file-backed mmap seam ------------------------------------------------
// File-backed mmap copies bytes from a host fd into the guest mapping, which
// requires host-fd access that lives only in the POSIX library. The portable
// SysMmap calls this hook for non-anonymous mappings; the POSIX TU registers a
// real loader at startup. In builds that link only the portable library
// (anonymous mappings / the heap use case) the hook stays null and file-backed
// mmap fails with EBADF. The loader returns std::nullopt on success or a Linux
// errno on failure; the caller owns unmapping on failure.
using MmapFileLoader = std::optional<int> (*)(ProcessId pid, AddressSpace& space, address_t dest, word_t length,
                                              word_t fd_arg, word_t offset) noexcept;
[[nodiscard]] MmapFileLoader mmap_file_loader() noexcept;
void set_mmap_file_loader(MmapFileLoader loader) noexcept;

} // namespace detail
} // namespace x86sim::linux_syscalls

#endif
