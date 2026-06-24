#include "x86sim-support/syscall-linux.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <filesystem>
#include <limits>
#include <netinet/in.h>
#include <new>
#include <span>
#include <string>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>
#include <utime.h>
#include <vector>

namespace x86sim::linux_syscalls {
namespace detail {

constexpr word_t syscall_read = 0;
constexpr word_t syscall_write = 1;
constexpr word_t syscall_open = 2;
constexpr word_t syscall_close = 3;
constexpr word_t syscall_stat = 4;
constexpr word_t syscall_fstat = 5;
constexpr word_t syscall_lstat = 6;
constexpr word_t syscall_lseek = 8;
constexpr word_t syscall_mmap = 9;
constexpr word_t syscall_munmap = 11;
constexpr word_t syscall_brk = 12;
constexpr word_t syscall_rt_sigaction = 13;
constexpr word_t syscall_rt_sigprocmask = 14;
constexpr word_t syscall_ioctl = 16;
constexpr word_t syscall_pread64 = 17;
constexpr word_t syscall_pwrite64 = 18;
constexpr word_t syscall_access = 21;
constexpr word_t syscall_pipe = 22;
constexpr word_t syscall_select = 23;
constexpr word_t syscall_mremap = 25;
constexpr word_t syscall_nanosleep = 35;
constexpr word_t syscall_clone = 56;
constexpr word_t syscall_fork = 57;
constexpr word_t syscall_vfork = 58;
constexpr word_t syscall_execve = 59;
constexpr word_t syscall_wait4 = 61;
constexpr word_t syscall_fcntl = 72;
constexpr word_t syscall_chdir = 80;
constexpr word_t syscall_fchdir = 81;
constexpr word_t syscall_rename = 82;
constexpr word_t syscall_mkdir = 83;
constexpr word_t syscall_rmdir = 84;
constexpr word_t syscall_unlink = 87;
constexpr word_t syscall_umask = 95;
constexpr word_t syscall_gettimeofday = 96;
constexpr word_t syscall_getpid = 39;
constexpr word_t syscall_socket = 41;
constexpr word_t syscall_connect = 42;
constexpr word_t syscall_exit = 60;
constexpr word_t syscall_uname = 63;
constexpr word_t syscall_getcwd = 79;
constexpr word_t syscall_readlink = 89;
constexpr word_t syscall_chmod = 90;
constexpr word_t syscall_fchmod = 91;
constexpr word_t syscall_getrusage = 77;
constexpr word_t syscall_getuid = 102;
constexpr word_t syscall_getgid = 104;
constexpr word_t syscall_geteuid = 107;
constexpr word_t syscall_getegid = 108;
constexpr word_t syscall_utime = 132;
constexpr word_t syscall_arch_prctl = 158;
constexpr word_t syscall_gettid = 186;
constexpr word_t syscall_futex = 202;
constexpr word_t syscall_getdents64 = 217;
constexpr word_t syscall_set_tid_address = 218;
constexpr word_t syscall_clock_gettime = 228;
constexpr word_t syscall_exit_group = 231;
constexpr word_t syscall_openat = 257;
constexpr word_t syscall_mkdirat = 258;
constexpr word_t syscall_newfstatat = 262;
constexpr word_t syscall_unlinkat = 263;
constexpr word_t syscall_renameat = 264;
constexpr word_t syscall_readlinkat = 267;
constexpr word_t syscall_faccessat = 269;
constexpr word_t syscall_set_robust_list = 273;
constexpr word_t syscall_prlimit64 = 302;
constexpr word_t syscall_rseq = 334;

constexpr int linux_esrch = 3;
constexpr int linux_eio = 5;
constexpr int linux_ebadf = 9;
constexpr int linux_echild = 10;
constexpr int linux_eagain = 11;
constexpr int linux_enomem = 12;
constexpr int linux_efault = 14;
constexpr int linux_enotty = 25;
constexpr int linux_einval = 22;
constexpr int linux_emfile = 24;
constexpr int linux_eafnosupport = 97;
constexpr int linux_enametoolong = 36;
constexpr int linux_enosys = 38;
constexpr int linux_etimedout = 110;

constexpr word_t linux_o_accmode = 00000003;
constexpr word_t linux_o_rdonly = 00000000;
constexpr word_t linux_o_wronly = 00000001;
constexpr word_t linux_o_rdwr = 00000002;
constexpr word_t linux_o_creat = 00000100;
constexpr word_t linux_o_excl = 00000200;
constexpr word_t linux_o_noctty = 00000400;
constexpr word_t linux_o_trunc = 00001000;
constexpr word_t linux_o_append = 00002000;
constexpr word_t linux_o_nonblock = 00004000;
constexpr word_t linux_o_dsync = 00010000;
constexpr word_t linux_o_largefile = 00100000;
constexpr word_t linux_o_directory = 00200000;
constexpr word_t linux_o_nofollow = 00400000;
constexpr word_t linux_o_cloexec = 02000000;
constexpr word_t linux_o_sync = 04010000;

constexpr word_t linux_at_fdcwd = static_cast<word_t>(-100);
constexpr word_t linux_at_symlink_nofollow = 0x100;
constexpr word_t linux_at_removedir = 0x200;
constexpr word_t linux_at_empty_path = 0x1000;

constexpr word_t linux_f_dupfd = 0;
constexpr word_t linux_f_getfd = 1;
constexpr word_t linux_f_setfd = 2;
constexpr word_t linux_f_getfl = 3;
constexpr word_t linux_f_setfl = 4;
constexpr word_t linux_f_getlk = 5;
constexpr word_t linux_f_setlk = 6;
constexpr word_t linux_f_setlkw = 7;
constexpr word_t linux_f_dupfd_cloexec = 1030;

constexpr word_t linux_fd_cloexec = 1;
constexpr word_t linux_f_rdlck = 0;
constexpr word_t linux_f_wrlck = 1;
constexpr word_t linux_f_unlck = 2;

constexpr word_t linux_prot_read = 0x1;
constexpr word_t linux_prot_write = 0x2;
constexpr word_t linux_prot_exec = 0x4;

constexpr word_t linux_map_shared = 0x01;
constexpr word_t linux_map_private = 0x02;
constexpr word_t linux_map_fixed = 0x10;
constexpr word_t linux_map_anonymous = 0x20;
constexpr word_t linux_map_denied_write = 0x800;
constexpr word_t linux_map_executable = 0x1000;
constexpr word_t linux_map_noreserve = 0x4000;
constexpr word_t linux_map_populate = 0x8000;
constexpr word_t linux_map_stack = 0x20000;

constexpr word_t linux_mremap_maymove = 1;
constexpr word_t linux_mremap_fixed = 2;
constexpr word_t linux_mremap_dontunmap = 4;

constexpr word_t linux_af_unix = 1;
constexpr word_t linux_af_inet = 2;
constexpr word_t linux_af_inet6 = 10;

constexpr word_t linux_sock_type_mask = 0xf;
constexpr word_t linux_sock_stream = 1;
constexpr word_t linux_sock_dgram = 2;
constexpr word_t linux_sock_raw = 3;
constexpr word_t linux_sock_seqpacket = 5;
constexpr word_t linux_sock_nonblock = 00004000;
constexpr word_t linux_sock_cloexec = 02000000;

constexpr word_t linux_arch_set_gs = 0x1001;
constexpr word_t linux_arch_set_fs = 0x1002;
constexpr word_t linux_arch_get_fs = 0x1003;
constexpr word_t linux_arch_get_gs = 0x1004;

constexpr word_t linux_futex_wait = 0;
constexpr word_t linux_futex_wake = 1;
constexpr word_t linux_futex_wait_bitset = 9;
constexpr word_t linux_futex_wake_bitset = 10;
constexpr word_t linux_futex_private_flag = 128;
constexpr word_t linux_futex_clock_realtime = 256;
constexpr word_t linux_futex_cmd_mask = ~(linux_futex_private_flag | linux_futex_clock_realtime);

constexpr word_t synthetic_pid = 1;
constexpr word_t synthetic_uid = 1000;
constexpr word_t synthetic_gid = 1000;

constexpr std::size_t max_path_length = 4096;
constexpr std::size_t io_chunk_size = 64 * 1024;

[[nodiscard]] word_t syscall_number(const RegisterFile& context) noexcept {
  return context[Register::rax];
}

[[nodiscard]] word_t syscall_arg(const RegisterFile& context, std::size_t index) noexcept {
  static constexpr Register registers[] = {Register::rdi, Register::rsi, Register::rdx,
                                           Register::r10, Register::r8,  Register::r9};
  return context[registers[index]];
}

[[nodiscard]] bool handles(const RegisterFile& context, SyscallKind kind, word_t number) noexcept {
  return kind == SyscallKind::syscall64 && syscall_number(context) == number;
}

[[nodiscard]] SyscallResult continue_result() {
  return {.reason = StopReason::host_request, .continue_execution = true, .message = {}};
}

[[nodiscard]] SyscallResult return_value(RegisterFile& context, std::int64_t value) {
  context[Register::rax] = static_cast<word_t>(value);
  return continue_result();
}

[[nodiscard]] SyscallResult return_error(RegisterFile& context, int error) {
  return return_value(context, -static_cast<std::int64_t>(error));
}

SyscallResult unsupported_syscall(const RegisterFile& context, SyscallKind kind) {
  return {.reason = StopReason::unsupported_syscall,
          .continue_execution = false,
          .message = "unsupported Linux syscall " + std::to_string(syscall_number(context)) + " via " +
                     (kind == SyscallKind::syscall64 ? "syscall64"
                      : kind == SyscallKind::int80   ? "int80"
                                                     : "sysenter")};
}

[[nodiscard]] int memory_error_to_linux(MemoryError error) noexcept {
  switch (error) {
  case MemoryError::unaligned_address:
  case MemoryError::zero_size:
    return linux_einval;
  case MemoryError::unmapped_address:
    return linux_efault;
  case MemoryError::out_of_memory:
  case MemoryError::mapping_failed:
    return linux_enomem;
  }
  return linux_eio;
}

[[nodiscard]] int host_errno_to_linux(int error) noexcept {
  switch (error) {
#ifdef EPERM
  case EPERM:
    return 1;
#endif
#ifdef ENOENT
  case ENOENT:
    return 2;
#endif
#ifdef ESRCH
  case ESRCH:
    return 3;
#endif
#ifdef EINTR
  case EINTR:
    return 4;
#endif
#ifdef EIO
  case EIO:
    return linux_eio;
#endif
#ifdef ENXIO
  case ENXIO:
    return 6;
#endif
#ifdef E2BIG
  case E2BIG:
    return 7;
#endif
#ifdef ENOEXEC
  case ENOEXEC:
    return 8;
#endif
#ifdef EBADF
  case EBADF:
    return linux_ebadf;
#endif
#ifdef ECHILD
  case ECHILD:
    return 10;
#endif
#ifdef EAGAIN
  case EAGAIN:
    return 11;
#endif
#ifdef ENOMEM
  case ENOMEM:
    return linux_enomem;
#endif
#ifdef EACCES
  case EACCES:
    return 13;
#endif
#ifdef EFAULT
  case EFAULT:
    return linux_efault;
#endif
#ifdef EBUSY
  case EBUSY:
    return 16;
#endif
#ifdef EEXIST
  case EEXIST:
    return 17;
#endif
#ifdef EXDEV
  case EXDEV:
    return 18;
#endif
#ifdef ENODEV
  case ENODEV:
    return 19;
#endif
#ifdef ENOTDIR
  case ENOTDIR:
    return 20;
#endif
#ifdef EISDIR
  case EISDIR:
    return 21;
#endif
#ifdef EINVAL
  case EINVAL:
    return linux_einval;
#endif
#ifdef ENFILE
  case ENFILE:
    return 23;
#endif
#ifdef EMFILE
  case EMFILE:
    return 24;
#endif
#ifdef ENOTTY
  case ENOTTY:
    return 25;
#endif
#ifdef ETXTBSY
  case ETXTBSY:
    return 26;
#endif
#ifdef EFBIG
  case EFBIG:
    return 27;
#endif
#ifdef ENOSPC
  case ENOSPC:
    return 28;
#endif
#ifdef ESPIPE
  case ESPIPE:
    return 29;
#endif
#ifdef EROFS
  case EROFS:
    return 30;
#endif
#ifdef EMLINK
  case EMLINK:
    return 31;
#endif
#ifdef EPIPE
  case EPIPE:
    return 32;
#endif
#ifdef EDOM
  case EDOM:
    return 33;
#endif
#ifdef ERANGE
  case ERANGE:
    return 34;
#endif
#ifdef ENOTEMPTY
  case ENOTEMPTY:
    return 39;
#endif
#ifdef ELOOP
  case ELOOP:
    return 40;
#endif
#ifdef ENAMETOOLONG
  case ENAMETOOLONG:
    return linux_enametoolong;
#endif
#ifdef ENOSYS
  case ENOSYS:
    return linux_enosys;
#endif
  default:
    return linux_eio;
  }
}

[[nodiscard]] std::optional<int> checked_fd(word_t raw_fd) noexcept {
  if (raw_fd > static_cast<word_t>(std::numeric_limits<int>::max()))
    return std::nullopt;
  return static_cast<int>(raw_fd);
}

[[nodiscard]] bool fits_host_transfer(word_t count) noexcept {
  return count <= static_cast<word_t>(std::numeric_limits<ssize_t>::max());
}

[[nodiscard]] bool range_overflows(address_t start, word_t length) noexcept {
  return length != 0 && start > std::numeric_limits<address_t>::max() - (length - 1);
}

[[nodiscard]] bool is_linux_at_fdcwd(word_t raw_fd) noexcept {
  return static_cast<std::int32_t>(raw_fd & 0xffffffffu) == -100;
}

[[nodiscard]] std::optional<int> dirfd_from_linux(word_t raw_fd) noexcept {
  if (is_linux_at_fdcwd(raw_fd))
    return AT_FDCWD;
  return checked_fd(raw_fd);
}

[[nodiscard]] std::optional<std::uint64_t> page_align_up(std::uint64_t value) noexcept {
  constexpr std::uint64_t mask = Machine::kPageSize - 1;
  if (value > std::numeric_limits<std::uint64_t>::max() - mask)
    return std::nullopt;
  return (value + mask) & ~mask;
}

[[nodiscard]] bool page_aligned(std::uint64_t value) noexcept {
  return (value & (Machine::kPageSize - 1)) == 0;
}

[[nodiscard]] bool range_is_unmapped(Machine& machine, address_t start, word_t length) noexcept {
  if (length == 0)
    return true;
  if (range_overflows(start, length))
    return false;

  std::array<std::byte, 1> byte{};
  const address_t end = start + length;
  for (address_t page = start; page < end;) {
    if (machine.read_memory_into(page, byte))
      return false;

    const address_t next_page = (page + Machine::kPageSize) & ~(Machine::kPageSize - 1);
    if (next_page <= page)
      return false;
    page = next_page;
  }
  return true;
}

[[nodiscard]] std::optional<address_t> find_unmapped_range(Machine& machine, address_t& cursor,
                                                           word_t length) noexcept {
  auto candidate = page_align_up(cursor);
  if (!candidate)
    return std::nullopt;

  while (*candidate <= std::numeric_limits<address_t>::max() - length) {
    if (range_is_unmapped(machine, *candidate, length)) {
      cursor = *candidate + length;
      return *candidate;
    }
    if (*candidate > std::numeric_limits<address_t>::max() - Machine::kPageSize)
      return std::nullopt;
    *candidate += Machine::kPageSize;
  }
  return std::nullopt;
}

[[nodiscard]] std::optional<Protection> protection_from_linux(word_t prot) noexcept {
  if ((prot & ~(linux_prot_read | linux_prot_write | linux_prot_exec)) != 0)
    return std::nullopt;

  Protection result = Protection::none;
  if ((prot & linux_prot_read) != 0)
    result = result | Protection::read;
  if ((prot & linux_prot_write) != 0)
    result = result | Protection::write;
  if ((prot & linux_prot_exec) != 0)
    result = result | Protection::execute;
  return result;
}

[[nodiscard]] std::optional<int> write_word(Machine& machine, address_t address, word_t value) noexcept {
  std::array<std::byte, sizeof(word_t)> bytes{};
  for (std::size_t i = 0; i < bytes.size(); ++i)
    bytes[i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);

  auto written = machine.write_memory(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

[[nodiscard]] std::optional<int> write_int32(Machine& machine, address_t address, std::int32_t value) noexcept {
  std::array<std::byte, sizeof(value)> bytes{};
  const auto raw = static_cast<std::uint32_t>(value);
  for (std::size_t i = 0; i < bytes.size(); ++i)
    bytes[i] = static_cast<std::byte>((raw >> (i * 8)) & 0xff);

  auto written = machine.write_memory(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

template<std::size_t Size>
[[nodiscard]] std::optional<int> read_guest_memory(Machine& machine, address_t address,
                                                   std::array<std::byte, Size>& out) noexcept {
  auto read = machine.read_memory_into(address, out);
  if (!read)
    return memory_error_to_linux(read.error());
  return std::nullopt;
}

template<std::size_t Size>
void write_le(std::array<std::byte, Size>& bytes, std::size_t offset, std::uint64_t value, std::size_t width) noexcept {
  for (std::size_t i = 0; i < width; ++i)
    bytes[offset + i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);
}

[[nodiscard]] std::uint64_t read_le(std::span<const std::byte> bytes, std::size_t offset, std::size_t width) noexcept {
  std::uint64_t value = 0;
  for (std::size_t i = 0; i < width; ++i)
    value |= static_cast<std::uint64_t>(std::to_integer<unsigned char>(bytes[offset + i])) << (i * 8);
  return value;
}

[[nodiscard]] SysPrlimit64::Limit read_rlimit64(std::span<const std::byte> bytes) noexcept {
  return {.current = read_le(bytes, 0, 8), .maximum = read_le(bytes, 8, 8)};
}

[[nodiscard]] std::optional<int> write_rlimit64(Machine& machine, address_t address,
                                                SysPrlimit64::Limit limit) noexcept {
  std::array<std::byte, 16> bytes{};
  write_le(bytes, 0, limit.current, 8);
  write_le(bytes, 8, limit.maximum, 8);

  auto written = machine.write_memory(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

[[nodiscard]] std::optional<int> write_flock64(Machine& machine, address_t address,
                                               const struct flock& host_lock) noexcept {
  std::array<std::byte, 32> bytes{};
  word_t type = linux_f_unlck;
  switch (host_lock.l_type) {
  case F_RDLCK:
    type = linux_f_rdlck;
    break;
  case F_WRLCK:
    type = linux_f_wrlck;
    break;
  case F_UNLCK:
    type = linux_f_unlck;
    break;
  default:
    return linux_einval;
  }

  write_le(bytes, 0, type, 2);
  write_le(bytes, 2, static_cast<std::uint64_t>(host_lock.l_whence), 2);
  write_le(bytes, 8, static_cast<std::uint64_t>(host_lock.l_start), 8);
  write_le(bytes, 16, static_cast<std::uint64_t>(host_lock.l_len), 8);
  write_le(bytes, 24, static_cast<std::uint64_t>(host_lock.l_pid), 4);

  auto written = machine.write_memory(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

[[nodiscard]] std::optional<int> write_stat(Machine& machine, address_t address,
                                            const struct stat& host_stat) noexcept {
  std::array<std::byte, 144> bytes{};
  write_le(bytes, 0, static_cast<std::uint64_t>(host_stat.st_dev), 8);
  write_le(bytes, 8, static_cast<std::uint64_t>(host_stat.st_ino), 8);
  write_le(bytes, 16, static_cast<std::uint64_t>(host_stat.st_nlink), 8);
  write_le(bytes, 24, static_cast<std::uint64_t>(host_stat.st_mode), 4);
  write_le(bytes, 28, static_cast<std::uint64_t>(host_stat.st_uid), 4);
  write_le(bytes, 32, static_cast<std::uint64_t>(host_stat.st_gid), 4);
  write_le(bytes, 40, static_cast<std::uint64_t>(host_stat.st_rdev), 8);
  write_le(bytes, 48, static_cast<std::uint64_t>(host_stat.st_size), 8);
  write_le(bytes, 56, static_cast<std::uint64_t>(host_stat.st_blksize), 8);
  write_le(bytes, 64, static_cast<std::uint64_t>(host_stat.st_blocks), 8);
#if defined(__APPLE__)
  write_le(bytes, 72, static_cast<std::uint64_t>(host_stat.st_atimespec.tv_sec), 8);
  write_le(bytes, 80, static_cast<std::uint64_t>(host_stat.st_atimespec.tv_nsec), 8);
  write_le(bytes, 88, static_cast<std::uint64_t>(host_stat.st_mtimespec.tv_sec), 8);
  write_le(bytes, 96, static_cast<std::uint64_t>(host_stat.st_mtimespec.tv_nsec), 8);
  write_le(bytes, 104, static_cast<std::uint64_t>(host_stat.st_ctimespec.tv_sec), 8);
  write_le(bytes, 112, static_cast<std::uint64_t>(host_stat.st_ctimespec.tv_nsec), 8);
#else
  write_le(bytes, 72, static_cast<std::uint64_t>(host_stat.st_atim.tv_sec), 8);
  write_le(bytes, 80, static_cast<std::uint64_t>(host_stat.st_atim.tv_nsec), 8);
  write_le(bytes, 88, static_cast<std::uint64_t>(host_stat.st_mtim.tv_sec), 8);
  write_le(bytes, 96, static_cast<std::uint64_t>(host_stat.st_mtim.tv_nsec), 8);
  write_le(bytes, 104, static_cast<std::uint64_t>(host_stat.st_ctim.tv_sec), 8);
  write_le(bytes, 112, static_cast<std::uint64_t>(host_stat.st_ctim.tv_nsec), 8);
#endif

  auto written = machine.write_memory(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

template<std::size_t Size>
void copy_c_string(std::array<std::byte, Size>& bytes, std::size_t offset, const char* value) noexcept {
  const std::size_t length = std::min<std::size_t>(std::strlen(value), 64);
  for (std::size_t i = 0; i < length; ++i)
    bytes[offset + i] = static_cast<std::byte>(value[i]);
}

[[nodiscard]] std::optional<int> translate_open_flags(word_t flags) noexcept {
  constexpr word_t supported_flags = linux_o_accmode | linux_o_creat | linux_o_excl | linux_o_noctty | linux_o_trunc |
                                     linux_o_append | linux_o_nonblock | linux_o_dsync | linux_o_largefile |
                                     linux_o_directory | linux_o_nofollow | linux_o_cloexec | linux_o_sync;
  if ((flags & ~supported_flags) != 0)
    return std::nullopt;

  int host_flags = 0;
  switch (flags & linux_o_accmode) {
  case linux_o_rdonly:
    host_flags |= O_RDONLY;
    break;
  case linux_o_wronly:
    host_flags |= O_WRONLY;
    break;
  case linux_o_rdwr:
    host_flags |= O_RDWR;
    break;
  default:
    return std::nullopt;
  }

  if ((flags & linux_o_creat) != 0)
    host_flags |= O_CREAT;
  if ((flags & linux_o_excl) != 0)
    host_flags |= O_EXCL;
#ifdef O_NOCTTY
  if ((flags & linux_o_noctty) != 0)
    host_flags |= O_NOCTTY;
#endif
  if ((flags & linux_o_trunc) != 0)
    host_flags |= O_TRUNC;
  if ((flags & linux_o_append) != 0)
    host_flags |= O_APPEND;
  if ((flags & linux_o_nonblock) != 0)
    host_flags |= O_NONBLOCK;
#ifdef O_DSYNC
  if ((flags & linux_o_dsync) != 0)
    host_flags |= O_DSYNC;
#endif
#ifdef O_DIRECTORY
  if ((flags & linux_o_directory) != 0)
    host_flags |= O_DIRECTORY;
#endif
#ifdef O_NOFOLLOW
  if ((flags & linux_o_nofollow) != 0)
    host_flags |= O_NOFOLLOW;
#endif
#ifdef O_CLOEXEC
  if ((flags & linux_o_cloexec) != 0)
    host_flags |= O_CLOEXEC;
#endif
#ifdef O_SYNC
  if ((flags & linux_o_sync) != 0)
    host_flags |= O_SYNC;
#endif

  return host_flags;
}

[[nodiscard]] word_t translate_host_status_flags(int flags) noexcept {
  word_t result = 0;
  switch (flags & O_ACCMODE) {
  case O_WRONLY:
    result |= linux_o_wronly;
    break;
  case O_RDWR:
    result |= linux_o_rdwr;
    break;
  default:
    result |= linux_o_rdonly;
    break;
  }
  if ((flags & O_APPEND) != 0)
    result |= linux_o_append;
#ifdef O_NONBLOCK
  if ((flags & O_NONBLOCK) != 0)
    result |= linux_o_nonblock;
#endif
#ifdef O_DSYNC
  if ((flags & O_DSYNC) != 0)
    result |= linux_o_dsync;
#endif
#ifdef O_SYNC
  if ((flags & O_SYNC) != 0)
    result |= linux_o_sync;
#endif
  return result;
}

[[nodiscard]] std::optional<int> translate_setfl_flags(word_t flags) noexcept {
  constexpr word_t mutable_flags = linux_o_append | linux_o_nonblock | linux_o_dsync | linux_o_sync;
  constexpr word_t supported_flags = mutable_flags | linux_o_accmode | linux_o_largefile;
  if ((flags & ~supported_flags) != 0)
    return std::nullopt;
  flags &= mutable_flags;

  int host_flags = 0;
  if ((flags & linux_o_append) != 0)
    host_flags |= O_APPEND;
#ifdef O_NONBLOCK
  if ((flags & linux_o_nonblock) != 0)
    host_flags |= O_NONBLOCK;
#endif
#ifdef O_DSYNC
  if ((flags & linux_o_dsync) != 0)
    host_flags |= O_DSYNC;
#endif
#ifdef O_SYNC
  if ((flags & linux_o_sync) != 0)
    host_flags |= O_SYNC;
#endif
  return host_flags;
}

[[nodiscard]] std::optional<struct flock> read_flock64(Machine& machine, address_t address, int& error) noexcept {
  error = 0;
  if (address == 0 || range_overflows(address, 32)) {
    error = linux_efault;
    return std::nullopt;
  }

  std::array<std::byte, 32> bytes{};
  if (auto read_error = read_guest_memory(machine, address, bytes)) {
    error = *read_error;
    return std::nullopt;
  }

  struct flock host_lock{};
  switch (read_le(bytes, 0, 2)) {
  case linux_f_rdlck:
    host_lock.l_type = F_RDLCK;
    break;
  case linux_f_wrlck:
    host_lock.l_type = F_WRLCK;
    break;
  case linux_f_unlck:
    host_lock.l_type = F_UNLCK;
    break;
  default:
    error = linux_einval;
    return std::nullopt;
  }

  switch (read_le(bytes, 2, 2)) {
  case 0:
    host_lock.l_whence = SEEK_SET;
    break;
  case 1:
    host_lock.l_whence = SEEK_CUR;
    break;
  case 2:
    host_lock.l_whence = SEEK_END;
    break;
  default:
    error = linux_einval;
    return std::nullopt;
  }

  host_lock.l_start = static_cast<off_t>(static_cast<std::int64_t>(read_le(bytes, 8, 8)));
  host_lock.l_len = static_cast<off_t>(static_cast<std::int64_t>(read_le(bytes, 16, 8)));
  host_lock.l_pid = static_cast<pid_t>(read_le(bytes, 24, 4));
  return host_lock;
}

[[nodiscard]] std::optional<int> translate_socket_domain(word_t domain) noexcept {
  switch (domain) {
  case linux_af_unix:
    return AF_UNIX;
  case linux_af_inet:
    return AF_INET;
#ifdef AF_INET6
  case linux_af_inet6:
    return AF_INET6;
#endif
  default:
    return std::nullopt;
  }
}

[[nodiscard]] std::optional<int> translate_socket_type(word_t type) noexcept {
  switch (type & linux_sock_type_mask) {
  case linux_sock_stream:
    return SOCK_STREAM;
  case linux_sock_dgram:
    return SOCK_DGRAM;
#ifdef SOCK_RAW
  case linux_sock_raw:
    return SOCK_RAW;
#endif
#ifdef SOCK_SEQPACKET
  case linux_sock_seqpacket:
    return SOCK_SEQPACKET;
#endif
  default:
    return std::nullopt;
  }
}

[[nodiscard]] unsigned char linux_dirent_type(unsigned char host_type) noexcept {
  switch (host_type) {
#ifdef DT_FIFO
  case DT_FIFO:
    return 1;
#endif
#ifdef DT_CHR
  case DT_CHR:
    return 2;
#endif
#ifdef DT_DIR
  case DT_DIR:
    return 4;
#endif
#ifdef DT_BLK
  case DT_BLK:
    return 6;
#endif
#ifdef DT_REG
  case DT_REG:
    return 8;
#endif
#ifdef DT_LNK
  case DT_LNK:
    return 10;
#endif
#ifdef DT_SOCK
  case DT_SOCK:
    return 12;
#endif
  default:
    return 0;
  }
}

struct CStringResult {
  bool ok = false;
  int error = 0;
  std::string value;
};

[[nodiscard]] CStringResult read_c_string(Machine& machine, address_t address,
                                          std::size_t max_length = max_path_length) noexcept {
  if (address == 0)
    return {.ok = false, .error = linux_efault, .value = {}};

  try {
    std::string result;
    result.reserve(128);

    std::array<std::byte, Machine::kPageSize> buffer{};
    std::size_t total = 0;
    while (total < max_length) {
      if (range_overflows(address, total + 1))
        return {.ok = false, .error = linux_efault, .value = {}};

      const address_t current = address + total;
      const std::size_t page_remaining = Machine::kPageSize - (current & (Machine::kPageSize - 1));
      const std::size_t chunk_size = std::min(max_length - total, page_remaining);
      if (range_overflows(address, total + chunk_size))
        return {.ok = false, .error = linux_efault, .value = {}};

      auto read = machine.read_memory_into(current, std::span<std::byte>(buffer.data(), chunk_size));
      if (!read)
        return {.ok = false, .error = memory_error_to_linux(read.error()), .value = {}};

      for (std::size_t i = 0; i < chunk_size; ++i) {
        const char ch = static_cast<char>(std::to_integer<unsigned char>(buffer[i]));
        if (ch == '\0')
          return {.ok = true, .error = 0, .value = std::move(result)};

        result.push_back(ch);
      }
      total += chunk_size;
    }
  } catch (const std::bad_alloc&) {
    return {.ok = false, .error = linux_enomem, .value = {}};
  } catch (...) {
    return {.ok = false, .error = linux_eio, .value = {}};
  }

  return {.ok = false, .error = linux_enametoolong, .value = {}};
}

} // namespace detail

std::optional<SyscallResult> SysRead::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_read))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const address_t buffer_address = detail::syscall_arg(context, 1);
  const word_t count = detail::syscall_arg(context, 2);
  if (!detail::fits_host_transfer(count))
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, count))
    return detail::return_error(context, detail::linux_efault);

  try {
    std::vector<std::byte> buffer(static_cast<std::size_t>(count));
    const ssize_t bytes_read = read(*fd, buffer.data(), buffer.size());
    if (bytes_read < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));

    if (bytes_read == 0)
      return detail::return_value(context, 0);

    auto written = machine.write_memory(
        buffer_address, std::span<const std::byte>(buffer.data(), static_cast<std::size_t>(bytes_read)));
    if (!written)
      return detail::return_error(context, detail::memory_error_to_linux(written.error()));

    return detail::return_value(context, bytes_read);
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysReadlink::try_syscall(Machine& machine, RegisterFile& context,
                                                      SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_readlink) &&
      !detail::handles(context, kind, detail::syscall_readlinkat)) {
    return std::nullopt;
  }

  const bool is_readlinkat = detail::syscall_number(context) == detail::syscall_readlinkat;
  const address_t path_address = detail::syscall_arg(context, is_readlinkat ? 1 : 0);
  const address_t buffer_address = detail::syscall_arg(context, is_readlinkat ? 2 : 1);
  const word_t buffer_size = detail::syscall_arg(context, is_readlinkat ? 3 : 2);
  if (!detail::fits_host_transfer(buffer_size))
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, buffer_size))
    return detail::return_error(context, detail::linux_efault);

  auto path = detail::read_c_string(machine, path_address);
  if (!path.ok)
    return detail::return_error(context, path.error);

  try {
    std::vector<std::byte> buffer(static_cast<std::size_t>(buffer_size));
    ssize_t bytes_read = -1;
    if (is_readlinkat) {
      const word_t raw_fd = detail::syscall_arg(context, 0);
      auto fd = detail::dirfd_from_linux(raw_fd);
      if (!fd)
        return detail::return_error(context, detail::linux_ebadf);
      bytes_read = readlinkat(*fd, path.value.c_str(), reinterpret_cast<char*>(buffer.data()), buffer.size());
    } else {
      bytes_read = readlink(path.value.c_str(), reinterpret_cast<char*>(buffer.data()), buffer.size());
    }
    if (bytes_read < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));

    auto written = machine.write_memory(
        buffer_address, std::span<const std::byte>(buffer.data(), static_cast<std::size_t>(bytes_read)));
    if (!written)
      return detail::return_error(context, detail::memory_error_to_linux(written.error()));

    return detail::return_value(context, bytes_read);
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysWrite::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_write))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const address_t buffer_address = detail::syscall_arg(context, 1);
  const word_t count = detail::syscall_arg(context, 2);
  if (!detail::fits_host_transfer(count))
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, count))
    return detail::return_error(context, detail::linux_efault);

  try {
    std::array<std::byte, Machine::kPageSize> buffer{};
    word_t total = 0;
    while (total < count) {
      const auto address = buffer_address + total;
      const auto page_remaining = static_cast<word_t>(Machine::kPageSize - (address & (Machine::kPageSize - 1)));
      const auto chunk_limit = std::min<word_t>(count - total, static_cast<word_t>(detail::io_chunk_size));
      const auto chunk_size = static_cast<std::size_t>(std::min(chunk_limit, page_remaining));
      auto read = machine.read_memory_into(address, std::span<std::byte>(buffer.data(), chunk_size));
      if (!read)
        return total == 0 ? detail::return_error(context, detail::memory_error_to_linux(read.error()))
                          : detail::return_value(context, static_cast<std::int64_t>(total));

      const ssize_t written = write(*fd, buffer.data(), chunk_size);
      if (written < 0)
        return total == 0 ? detail::return_error(context, detail::host_errno_to_linux(errno))
                          : detail::return_value(context, static_cast<std::int64_t>(total));

      total += static_cast<word_t>(written);
      if (static_cast<std::size_t>(written) < chunk_size)
        return detail::return_value(context, static_cast<std::int64_t>(total));
    }

    return detail::return_value(context, static_cast<std::int64_t>(total));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysPreadPwrite::try_syscall(Machine& machine, RegisterFile& context,
                                                         SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_pread64) &&
      !detail::handles(context, kind, detail::syscall_pwrite64)) {
    return std::nullopt;
  }

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const bool is_pwrite = detail::syscall_number(context) == detail::syscall_pwrite64;
  const address_t buffer_address = detail::syscall_arg(context, 1);
  const word_t count = detail::syscall_arg(context, 2);
  const auto offset = static_cast<off_t>(detail::syscall_arg(context, 3));
  if (!detail::fits_host_transfer(count))
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, count))
    return detail::return_error(context, detail::linux_efault);

  try {
    if (!is_pwrite) {
      std::vector<std::byte> buffer(static_cast<std::size_t>(count));
      const ssize_t bytes_read = pread(*fd, buffer.data(), buffer.size(), offset);
      if (bytes_read < 0)
        return detail::return_error(context, detail::host_errno_to_linux(errno));

      if (bytes_read == 0)
        return detail::return_value(context, 0);

      auto written = machine.write_memory(
          buffer_address, std::span<const std::byte>(buffer.data(), static_cast<std::size_t>(bytes_read)));
      if (!written)
        return detail::return_error(context, detail::memory_error_to_linux(written.error()));
      return detail::return_value(context, bytes_read);
    }

    std::array<std::byte, Machine::kPageSize> buffer{};
    word_t total = 0;
    while (total < count) {
      const auto address = buffer_address + total;
      const auto page_remaining = static_cast<word_t>(Machine::kPageSize - (address & (Machine::kPageSize - 1)));
      const auto chunk_limit = std::min<word_t>(count - total, static_cast<word_t>(detail::io_chunk_size));
      const auto chunk_size = static_cast<std::size_t>(std::min(chunk_limit, page_remaining));
      auto read = machine.read_memory_into(address, std::span<std::byte>(buffer.data(), chunk_size));
      if (!read)
        return total == 0 ? detail::return_error(context, detail::memory_error_to_linux(read.error()))
                          : detail::return_value(context, static_cast<std::int64_t>(total));

      const ssize_t written = pwrite(*fd, buffer.data(), chunk_size, offset + static_cast<off_t>(total));
      if (written < 0)
        return total == 0 ? detail::return_error(context, detail::host_errno_to_linux(errno))
                          : detail::return_value(context, static_cast<std::int64_t>(total));

      total += static_cast<word_t>(written);
      if (static_cast<std::size_t>(written) < chunk_size)
        return detail::return_value(context, static_cast<std::int64_t>(total));
    }

    return detail::return_value(context, static_cast<std::int64_t>(total));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysOpen::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_open) && !detail::handles(context, kind, detail::syscall_openat))
    return std::nullopt;

  const bool is_openat = detail::syscall_number(context) == detail::syscall_openat;
  const address_t path_address = detail::syscall_arg(context, is_openat ? 1 : 0);
  auto path = detail::read_c_string(machine, path_address);
  if (!path.ok)
    return detail::return_error(context, path.error);

  auto flags = detail::translate_open_flags(detail::syscall_arg(context, is_openat ? 2 : 1));
  if (!flags)
    return detail::return_error(context, detail::linux_einval);

  const auto mode = static_cast<mode_t>(detail::syscall_arg(context, is_openat ? 3 : 2) & 07777);
  int fd = -1;
  if (is_openat) {
    const word_t raw_dirfd = detail::syscall_arg(context, 0);
    auto dirfd = detail::dirfd_from_linux(raw_dirfd);
    if (!dirfd)
      return detail::return_error(context, detail::linux_ebadf);
    fd = openat(*dirfd, path.value.c_str(), *flags, mode);
  } else {
    fd = open(path.value.c_str(), *flags, mode);
  }
  if (fd < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  return detail::return_value(context, fd);
}

std::optional<SyscallResult> SysClose::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_close))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  if (close(*fd) < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysPipe::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_pipe))
    return std::nullopt;

  const address_t pipefd_address = detail::syscall_arg(context, 0);
  if (pipefd_address == 0 || detail::range_overflows(pipefd_address, 8))
    return detail::return_error(context, detail::linux_efault);

  int pipefds[2] = {-1, -1};
  if (pipe(pipefds) < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  if (auto error = detail::write_int32(machine, pipefd_address, pipefds[0])) {
    close(pipefds[0]);
    close(pipefds[1]);
    return detail::return_error(context, *error);
  }
  if (auto error = detail::write_int32(machine, pipefd_address + 4, pipefds[1])) {
    close(pipefds[0]);
    close(pipefds[1]);
    return detail::return_error(context, *error);
  }

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysFstat::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_fstat))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  struct stat host_stat{};
  if (fstat(*fd, &host_stat) < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  if (auto error = detail::write_stat(machine, detail::syscall_arg(context, 1), host_stat))
    return detail::return_error(context, *error);
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysStat::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_stat) && !detail::handles(context, kind, detail::syscall_lstat) &&
      !detail::handles(context, kind, detail::syscall_newfstatat)) {
    return std::nullopt;
  }

  const word_t number = detail::syscall_number(context);
  const bool is_newfstatat = number == detail::syscall_newfstatat;
  const address_t path_address = detail::syscall_arg(context, is_newfstatat ? 1 : 0);
  const address_t stat_address = detail::syscall_arg(context, is_newfstatat ? 2 : 1);
  const word_t flags = is_newfstatat ? detail::syscall_arg(context, 3) : 0;
  constexpr word_t supported_flags = detail::linux_at_symlink_nofollow | detail::linux_at_empty_path;
  if ((flags & ~supported_flags) != 0)
    return detail::return_error(context, detail::linux_einval);

  auto path = detail::read_c_string(machine, path_address);
  if (!path.ok)
    return detail::return_error(context, path.error);

  struct stat host_stat{};
  int result = -1;
  if (is_newfstatat) {
    const word_t raw_dirfd = detail::syscall_arg(context, 0);
    auto dirfd = detail::dirfd_from_linux(raw_dirfd);
    if (!dirfd)
      return detail::return_error(context, detail::linux_ebadf);
    int host_flags = 0;
#ifdef AT_SYMLINK_NOFOLLOW
    if ((flags & detail::linux_at_symlink_nofollow) != 0)
      host_flags |= AT_SYMLINK_NOFOLLOW;
#endif
#ifdef AT_EMPTY_PATH
    if ((flags & detail::linux_at_empty_path) != 0)
      host_flags |= AT_EMPTY_PATH;
#endif
    result = fstatat(*dirfd, path.value.c_str(), &host_stat, host_flags);
  } else if (number == detail::syscall_lstat) {
    result = lstat(path.value.c_str(), &host_stat);
  } else {
    result = stat(path.value.c_str(), &host_stat);
  }

  if (result < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  if (auto error = detail::write_stat(machine, stat_address, host_stat))
    return detail::return_error(context, *error);
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysGetdents64::try_syscall(Machine& machine, RegisterFile& context,
                                                        SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_getdents64))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const address_t buffer_address = detail::syscall_arg(context, 1);
  const word_t count = detail::syscall_arg(context, 2);
  if (!detail::fits_host_transfer(count))
    return detail::return_error(context, detail::linux_einval);
  if (buffer_address == 0 || detail::range_overflows(buffer_address, count))
    return detail::return_error(context, detail::linux_efault);

  const int duplicated_fd = dup(*fd);
  if (duplicated_fd < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  DIR* directory = fdopendir(duplicated_fd);
  if (!directory) {
    const int saved_errno = errno;
    close(duplicated_fd);
    return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
  }

  try {
    std::vector<std::byte> output;
    output.reserve(static_cast<std::size_t>(std::min<word_t>(count, detail::io_chunk_size)));

    while (output.size() < count) {
      const long before = telldir(directory);
      errno = 0;
      dirent* entry = readdir(directory);
      if (!entry) {
        if (errno != 0) {
          const int saved_errno = errno;
          closedir(directory);
          return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
        }
        break;
      }

      const std::size_t name_length = std::strlen(entry->d_name);
      const std::size_t record_length = (19 + name_length + 1 + 7) & ~std::size_t{7};
      if (record_length > count - output.size()) {
        seekdir(directory, before);
        break;
      }

      const long next = telldir(directory);
      const std::size_t offset = output.size();
      output.resize(output.size() + record_length);
      std::fill(output.begin() + static_cast<std::ptrdiff_t>(offset), output.end(), std::byte{0});

      auto write_field = [&](std::size_t field_offset, std::uint64_t value, std::size_t width) noexcept {
        for (std::size_t i = 0; i < width; ++i)
          output[offset + field_offset + i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);
      };

      write_field(0, static_cast<std::uint64_t>(entry->d_ino), 8);
      write_field(8, static_cast<std::uint64_t>(next), 8);
      write_field(16, static_cast<std::uint64_t>(record_length), 2);
#ifdef DT_UNKNOWN
      output[offset + 18] = static_cast<std::byte>(detail::linux_dirent_type(entry->d_type));
#else
      output[offset + 18] = std::byte{0};
#endif
      std::memcpy(output.data() + offset + 19, entry->d_name, name_length);
    }

    if (output.empty()) {
      closedir(directory);
      return detail::return_value(context, 0);
    }

    auto written = machine.write_memory(buffer_address, output);
    if (!written) {
      closedir(directory);
      return detail::return_error(context, detail::memory_error_to_linux(written.error()));
    }

    const auto bytes_written = static_cast<std::int64_t>(output.size());
    closedir(directory);
    return detail::return_value(context, bytes_written);
  } catch (const std::bad_alloc&) {
    closedir(directory);
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    closedir(directory);
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysFileSystem::try_syscall(Machine& machine, RegisterFile& context,
                                                        SyscallKind kind) noexcept {
  const word_t number = detail::syscall_number(context);
  if (!detail::handles(context, kind, detail::syscall_access) &&
      !detail::handles(context, kind, detail::syscall_chdir) &&
      !detail::handles(context, kind, detail::syscall_fchdir) &&
      !detail::handles(context, kind, detail::syscall_rename) &&
      !detail::handles(context, kind, detail::syscall_mkdir) &&
      !detail::handles(context, kind, detail::syscall_rmdir) &&
      !detail::handles(context, kind, detail::syscall_unlink) &&
      !detail::handles(context, kind, detail::syscall_chmod) &&
      !detail::handles(context, kind, detail::syscall_fchmod) &&
      !detail::handles(context, kind, detail::syscall_umask) &&
      !detail::handles(context, kind, detail::syscall_utime) &&
      !detail::handles(context, kind, detail::syscall_mkdirat) &&
      !detail::handles(context, kind, detail::syscall_unlinkat) &&
      !detail::handles(context, kind, detail::syscall_renameat) &&
      !detail::handles(context, kind, detail::syscall_faccessat)) {
    return std::nullopt;
  }

  switch (number) {
  case detail::syscall_access: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    if (access(path.value.c_str(), static_cast<int>(detail::syscall_arg(context, 1))) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_chdir: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    if (chdir(path.value.c_str()) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_fchdir: {
    auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
    if (!fd)
      return detail::return_error(context, detail::linux_ebadf);

    if (fchdir(*fd) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_rename: {
    auto old_path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!old_path.ok)
      return detail::return_error(context, old_path.error);
    auto new_path = detail::read_c_string(machine, detail::syscall_arg(context, 1));
    if (!new_path.ok)
      return detail::return_error(context, new_path.error);

    if (rename(old_path.value.c_str(), new_path.value.c_str()) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_mkdir: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const auto mode = static_cast<mode_t>(detail::syscall_arg(context, 1) & 07777);
    if (mkdir(path.value.c_str(), mode) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_rmdir: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    if (rmdir(path.value.c_str()) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_unlink: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    if (unlink(path.value.c_str()) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_chmod: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const auto mode = static_cast<mode_t>(detail::syscall_arg(context, 1) & 07777);
    if (chmod(path.value.c_str(), mode) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_fchmod: {
    auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
    if (!fd)
      return detail::return_error(context, detail::linux_ebadf);

    const auto mode = static_cast<mode_t>(detail::syscall_arg(context, 1) & 07777);
    if (fchmod(*fd, mode) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_umask: {
    const mode_t old_mask = umask(static_cast<mode_t>(detail::syscall_arg(context, 0) & 0777));
    return detail::return_value(context, old_mask);
  }
  case detail::syscall_utime: {
    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 0));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const address_t times_address = detail::syscall_arg(context, 1);
    if (times_address == 0) {
      if (utime(path.value.c_str(), nullptr) < 0)
        return detail::return_error(context, detail::host_errno_to_linux(errno));
      return detail::return_value(context, 0);
    }

    std::array<std::byte, 16> bytes{};
    if (auto read_error = detail::read_guest_memory(machine, times_address, bytes))
      return detail::return_error(context, *read_error);

    utimbuf times{};
    times.actime = static_cast<time_t>(static_cast<std::int64_t>(detail::read_le(bytes, 0, 8)));
    times.modtime = static_cast<time_t>(static_cast<std::int64_t>(detail::read_le(bytes, 8, 8)));
    if (utime(path.value.c_str(), &times) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_mkdirat: {
    auto dirfd = detail::dirfd_from_linux(detail::syscall_arg(context, 0));
    if (!dirfd)
      return detail::return_error(context, detail::linux_ebadf);

    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 1));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const auto mode = static_cast<mode_t>(detail::syscall_arg(context, 2) & 07777);
    if (mkdirat(*dirfd, path.value.c_str(), mode) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_unlinkat: {
    auto dirfd = detail::dirfd_from_linux(detail::syscall_arg(context, 0));
    if (!dirfd)
      return detail::return_error(context, detail::linux_ebadf);

    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 1));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const word_t flags = detail::syscall_arg(context, 2);
    if ((flags & ~detail::linux_at_removedir) != 0)
      return detail::return_error(context, detail::linux_einval);

    if (unlinkat(*dirfd, path.value.c_str(), (flags & detail::linux_at_removedir) ? AT_REMOVEDIR : 0) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_renameat: {
    auto old_dirfd = detail::dirfd_from_linux(detail::syscall_arg(context, 0));
    auto new_dirfd = detail::dirfd_from_linux(detail::syscall_arg(context, 2));
    if (!old_dirfd || !new_dirfd)
      return detail::return_error(context, detail::linux_ebadf);

    auto old_path = detail::read_c_string(machine, detail::syscall_arg(context, 1));
    if (!old_path.ok)
      return detail::return_error(context, old_path.error);
    auto new_path = detail::read_c_string(machine, detail::syscall_arg(context, 3));
    if (!new_path.ok)
      return detail::return_error(context, new_path.error);

    if (renameat(*old_dirfd, old_path.value.c_str(), *new_dirfd, new_path.value.c_str()) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::syscall_faccessat: {
    auto dirfd = detail::dirfd_from_linux(detail::syscall_arg(context, 0));
    if (!dirfd)
      return detail::return_error(context, detail::linux_ebadf);

    auto path = detail::read_c_string(machine, detail::syscall_arg(context, 1));
    if (!path.ok)
      return detail::return_error(context, path.error);

    const int mode = static_cast<int>(detail::syscall_arg(context, 2));
    const int flags = static_cast<int>(detail::syscall_arg(context, 3));
    if (faccessat(*dirfd, path.value.c_str(), mode, flags) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  }

  return std::nullopt;
}

std::optional<SyscallResult> SysLseek::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_lseek))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const auto offset = static_cast<off_t>(detail::syscall_arg(context, 1));
  int whence = 0;
  switch (detail::syscall_arg(context, 2)) {
  case 0:
    whence = SEEK_SET;
    break;
  case 1:
    whence = SEEK_CUR;
    break;
  case 2:
    whence = SEEK_END;
    break;
  default:
    return detail::return_error(context, detail::linux_einval);
  }

  const off_t result = lseek(*fd, offset, whence);
  if (result < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));
  return detail::return_value(context, result);
}

std::optional<SyscallResult> SysIoctl::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_ioctl))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  return detail::return_error(context, detail::linux_enotty);
}

std::optional<SyscallResult> SysFcntl::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_fcntl))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const word_t command = detail::syscall_arg(context, 1);
  const word_t arg = detail::syscall_arg(context, 2);

  switch (command) {
  case detail::linux_f_dupfd:
  case detail::linux_f_dupfd_cloexec: {
    if (arg > static_cast<word_t>(std::numeric_limits<int>::max()))
      return detail::return_error(context, detail::linux_einval);
    const int duplicated = fcntl(*fd, F_DUPFD, static_cast<int>(arg));
    if (duplicated < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    if (command == detail::linux_f_dupfd_cloexec && fcntl(duplicated, F_SETFD, FD_CLOEXEC) < 0) {
      const int saved_errno = errno;
      close(duplicated);
      return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
    }
    return detail::return_value(context, duplicated);
  }
  case detail::linux_f_getfd: {
    const int flags = fcntl(*fd, F_GETFD, 0);
    if (flags < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, (flags & FD_CLOEXEC) != 0 ? detail::linux_fd_cloexec : 0);
  }
  case detail::linux_f_setfd: {
    const int flags = (arg & detail::linux_fd_cloexec) != 0 ? FD_CLOEXEC : 0;
    if (fcntl(*fd, F_SETFD, flags) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::linux_f_getfl: {
    const int flags = fcntl(*fd, F_GETFL, 0);
    if (flags < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, static_cast<std::int64_t>(detail::translate_host_status_flags(flags)));
  }
  case detail::linux_f_setfl: {
    const auto flags = detail::translate_setfl_flags(arg);
    if (!flags)
      return detail::return_error(context, detail::linux_einval);
    if (fcntl(*fd, F_SETFL, *flags) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }
  case detail::linux_f_getlk:
  case detail::linux_f_setlk:
  case detail::linux_f_setlkw: {
    int error = 0;
    auto host_lock = detail::read_flock64(machine, arg, error);
    if (!host_lock)
      return detail::return_error(context, error);

    int host_command = F_GETLK;
    if (command == detail::linux_f_setlk)
      host_command = F_SETLK;
    else if (command == detail::linux_f_setlkw)
      host_command = F_SETLKW;

    if (fcntl(*fd, host_command, &*host_lock) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));

    if (command == detail::linux_f_getlk) {
      if (auto write_error = detail::write_flock64(machine, arg, *host_lock))
        return detail::return_error(context, *write_error);
    }
    return detail::return_value(context, 0);
  }
  default:
    return detail::return_error(context, detail::linux_einval);
  }
}

std::optional<SyscallResult> SysSocket::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_socket))
    return std::nullopt;

  const word_t raw_domain = detail::syscall_arg(context, 0);
  const word_t raw_type = detail::syscall_arg(context, 1);
  const word_t raw_protocol = detail::syscall_arg(context, 2);
  constexpr word_t supported_flags = detail::linux_sock_nonblock | detail::linux_sock_cloexec;
  if ((raw_type & ~(detail::linux_sock_type_mask | supported_flags)) != 0)
    return detail::return_error(context, detail::linux_einval);

  const auto domain = detail::translate_socket_domain(raw_domain);
  if (!domain)
    return detail::return_error(context, detail::linux_eafnosupport);

  const auto type = detail::translate_socket_type(raw_type);
  if (!type)
    return detail::return_error(context, detail::linux_einval);

  if (raw_protocol > static_cast<word_t>(std::numeric_limits<int>::max()))
    return detail::return_error(context, detail::linux_einval);

  const int fd = socket(*domain, *type, static_cast<int>(raw_protocol));
  if (fd < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  if ((raw_type & detail::linux_sock_cloexec) != 0) {
    if (fcntl(fd, F_SETFD, FD_CLOEXEC) < 0) {
      const int saved_errno = errno;
      close(fd);
      return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
    }
  }

  if ((raw_type & detail::linux_sock_nonblock) != 0) {
    const int flags = fcntl(fd, F_GETFL, 0);
    if (flags < 0 || fcntl(fd, F_SETFL, flags | O_NONBLOCK) < 0) {
      const int saved_errno = errno;
      close(fd);
      return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
    }
  }

  return detail::return_value(context, fd);
}

std::optional<SyscallResult> SysConnect::try_syscall(Machine& machine, RegisterFile& context,
                                                     SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_connect))
    return std::nullopt;

  const auto fd = detail::checked_fd(detail::syscall_arg(context, 0));
  if (!fd)
    return detail::return_error(context, detail::linux_ebadf);

  const address_t address = detail::syscall_arg(context, 1);
  const word_t length = detail::syscall_arg(context, 2);
  if (address == 0 || length < 2 || length > 4096 || detail::range_overflows(address, length))
    return detail::return_error(context, detail::linux_efault);

  try {
    std::vector<std::byte> bytes(static_cast<std::size_t>(length));
    auto read = machine.read_memory_into(address, bytes);
    if (!read)
      return detail::return_error(context, detail::memory_error_to_linux(read.error()));

    const word_t family = detail::read_le(bytes, 0, 2);
    int result = -1;

    switch (family) {
    case detail::linux_af_unix: {
      if (length < 3)
        return detail::return_error(context, detail::linux_einval);

      sockaddr_un host_address{};
#if defined(__APPLE__)
      host_address.sun_len = sizeof(sockaddr_un);
#endif
      host_address.sun_family = AF_UNIX;

      const std::size_t path_length = std::min<std::size_t>(bytes.size() - 2, sizeof(host_address.sun_path) - 1);
      for (std::size_t i = 0; i < path_length; ++i)
        host_address.sun_path[i] = static_cast<char>(std::to_integer<unsigned char>(bytes[2 + i]));

      if (host_address.sun_path[0] == '\0')
        return detail::return_error(context, detail::linux_eafnosupport);

      result = connect(*fd, reinterpret_cast<sockaddr*>(&host_address), sizeof(host_address));
      break;
    }
    case detail::linux_af_inet: {
      if (length < 16)
        return detail::return_error(context, detail::linux_einval);

      sockaddr_in host_address{};
#if defined(__APPLE__)
      host_address.sin_len = sizeof(sockaddr_in);
#endif
      host_address.sin_family = AF_INET;
      host_address.sin_port = static_cast<in_port_t>(detail::read_le(bytes, 2, 2));
      host_address.sin_addr.s_addr = static_cast<in_addr_t>(detail::read_le(bytes, 4, 4));
      result = connect(*fd, reinterpret_cast<sockaddr*>(&host_address), sizeof(host_address));
      break;
    }
#ifdef AF_INET6
    case detail::linux_af_inet6: {
      if (length < 28)
        return detail::return_error(context, detail::linux_einval);

      sockaddr_in6 host_address{};
#if defined(__APPLE__)
      host_address.sin6_len = sizeof(sockaddr_in6);
#endif
      host_address.sin6_family = AF_INET6;
      host_address.sin6_port = static_cast<in_port_t>(detail::read_le(bytes, 2, 2));
      host_address.sin6_flowinfo = static_cast<std::uint32_t>(detail::read_le(bytes, 4, 4));
      for (std::size_t i = 0; i < 16; ++i)
        host_address.sin6_addr.s6_addr[i] = std::to_integer<unsigned char>(bytes[8 + i]);
      host_address.sin6_scope_id = static_cast<std::uint32_t>(detail::read_le(bytes, 24, 4));
      result = connect(*fd, reinterpret_cast<sockaddr*>(&host_address), sizeof(host_address));
      break;
    }
#endif
    default:
      return detail::return_error(context, detail::linux_eafnosupport);
    }

    if (result < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysSelect::try_syscall(Machine& machine, RegisterFile& context,
                                                    SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_select))
    return std::nullopt;

  const word_t nfds = detail::syscall_arg(context, 0);
  const address_t readfds = detail::syscall_arg(context, 1);
  const address_t writefds = detail::syscall_arg(context, 2);
  const address_t exceptfds = detail::syscall_arg(context, 3);
  const address_t timeout_address = detail::syscall_arg(context, 4);

  if (nfds > static_cast<word_t>(std::numeric_limits<int>::max()))
    return detail::return_error(context, detail::linux_einval);

  const word_t fdset_bytes = (nfds + 7) / 8;
  if ((readfds != 0 && detail::range_overflows(readfds, fdset_bytes)) ||
      (writefds != 0 && detail::range_overflows(writefds, fdset_bytes)) ||
      (exceptfds != 0 && detail::range_overflows(exceptfds, fdset_bytes))) {
    return detail::return_error(context, detail::linux_efault);
  }

  if (timeout_address != 0) {
    std::array<std::byte, 16> bytes{};
    if (auto read_error = detail::read_guest_memory(machine, timeout_address, bytes))
      return detail::return_error(context, *read_error);

    timespec requested{};
    requested.tv_sec = static_cast<time_t>(static_cast<std::int64_t>(detail::read_le(bytes, 0, 8)));
    requested.tv_nsec = static_cast<long>(detail::read_le(bytes, 8, 8)) * 1000L;
    if (requested.tv_sec < 0 || requested.tv_nsec < 0 || requested.tv_nsec >= 1000000000L)
      return detail::return_error(context, detail::linux_einval);
    if (requested.tv_sec != 0 || requested.tv_nsec != 0)
      nanosleep(&requested, nullptr);
  }

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysProcess::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_clone) && !detail::handles(context, kind, detail::syscall_fork) &&
      !detail::handles(context, kind, detail::syscall_vfork) &&
      !detail::handles(context, kind, detail::syscall_execve) &&
      !detail::handles(context, kind, detail::syscall_wait4)) {
    return std::nullopt;
  }

  if (detail::syscall_number(context) == detail::syscall_wait4)
    return detail::return_error(context, detail::linux_echild);
  return detail::return_error(context, detail::linux_enosys);
}

std::optional<SyscallResult> SysFutex::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_futex))
    return std::nullopt;

  const address_t uaddr = detail::syscall_arg(context, 0);
  const word_t op = detail::syscall_arg(context, 1);
  const std::uint32_t expected = static_cast<std::uint32_t>(detail::syscall_arg(context, 2));
  const address_t timeout = detail::syscall_arg(context, 3);
  const word_t command = op & detail::linux_futex_cmd_mask;

  switch (command) {
  case detail::linux_futex_wait:
  case detail::linux_futex_wait_bitset: {
    if (uaddr == 0 || detail::range_overflows(uaddr, 4))
      return detail::return_error(context, detail::linux_efault);

    std::array<std::byte, 4> bytes{};
    if (auto read_error = detail::read_guest_memory(machine, uaddr, bytes))
      return detail::return_error(context, *read_error);

    const auto current = static_cast<std::uint32_t>(detail::read_le(bytes, 0, 4));
    if (current != expected)
      return detail::return_error(context, detail::linux_eagain);

    if (timeout != 0)
      return detail::return_error(context, detail::linux_etimedout);

    return detail::return_value(context, 0);
  }
  case detail::linux_futex_wake:
  case detail::linux_futex_wake_bitset:
    return detail::return_value(context, 0);
  default:
    return detail::return_error(context, detail::linux_enosys);
  }
}

std::optional<SyscallResult> SysSignals::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_rt_sigaction) &&
      !detail::handles(context, kind, detail::syscall_rt_sigprocmask)) {
    return std::nullopt;
  }

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysBrk::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_brk))
    return std::nullopt;

  const address_t requested_break = detail::syscall_arg(context, 0);
  if (requested_break == 0)
    return detail::return_value(context, static_cast<std::int64_t>(current_break));
  if (requested_break < minimum_break)
    return detail::return_value(context, static_cast<std::int64_t>(current_break));

  auto requested_end = detail::page_align_up(requested_break);
  if (!requested_end)
    return detail::return_value(context, static_cast<std::int64_t>(current_break));

  if (*requested_end > mapped_end) {
    auto mapped = machine.map(mapped_end, *requested_end - mapped_end, Protection::read | Protection::write);
    if (!mapped)
      return detail::return_value(context, static_cast<std::int64_t>(current_break));
    mapped_end = *requested_end;
  }

  current_break = requested_break;
  return detail::return_value(context, static_cast<std::int64_t>(current_break));
}

std::optional<SyscallResult> SysArchPrctl::try_syscall(Machine& machine, RegisterFile& context,
                                                       SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_arch_prctl))
    return std::nullopt;

  const word_t code = detail::syscall_arg(context, 0);
  const address_t address = detail::syscall_arg(context, 1);

  switch (code) {
  case detail::linux_arch_set_fs:
    machine.set_segment_base(context, SegmentRegister::fs, address);
    return detail::return_value(context, 0);
  case detail::linux_arch_set_gs:
    machine.set_segment_base(context, SegmentRegister::gs, address);
    return detail::return_value(context, 0);
  case detail::linux_arch_get_fs:
    if (auto error = detail::write_word(machine, address, machine.segment_base(context, SegmentRegister::fs)))
      return detail::return_error(context, *error);
    return detail::return_value(context, 0);
  case detail::linux_arch_get_gs:
    if (auto error = detail::write_word(machine, address, machine.segment_base(context, SegmentRegister::gs)))
      return detail::return_error(context, *error);
    return detail::return_value(context, 0);
  default:
    return detail::return_error(context, detail::linux_einval);
  }
}

std::optional<SyscallResult> SysGetIdentity::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_getpid) &&
      !detail::handles(context, kind, detail::syscall_gettid) &&
      !detail::handles(context, kind, detail::syscall_getuid) &&
      !detail::handles(context, kind, detail::syscall_geteuid) &&
      !detail::handles(context, kind, detail::syscall_getgid) &&
      !detail::handles(context, kind, detail::syscall_getegid)) {
    return std::nullopt;
  }

  switch (detail::syscall_number(context)) {
  case detail::syscall_getpid:
  case detail::syscall_gettid:
    return detail::return_value(context, detail::synthetic_pid);
  case detail::syscall_getuid:
  case detail::syscall_geteuid:
    return detail::return_value(context, detail::synthetic_uid);
  case detail::syscall_getgid:
  case detail::syscall_getegid:
    return detail::return_value(context, detail::synthetic_gid);
  default:
    return std::nullopt;
  }
}

std::optional<SyscallResult> SysUname::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_uname))
    return std::nullopt;

  constexpr std::size_t field_size = 65;
  std::array<std::byte, field_size * 6> utsname{};
  detail::copy_c_string(utsname, field_size * 0, "Linux");
  detail::copy_c_string(utsname, field_size * 1, "x86sim");
  detail::copy_c_string(utsname, field_size * 2, "6.0.0");
  detail::copy_c_string(utsname, field_size * 3, "#1");
  detail::copy_c_string(utsname, field_size * 4, "x86_64");
  detail::copy_c_string(utsname, field_size * 5, "localdomain");

  auto written = machine.write_memory(detail::syscall_arg(context, 0), utsname);
  if (!written)
    return detail::return_error(context, detail::memory_error_to_linux(written.error()));
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysGetcwd::try_syscall(Machine& machine, RegisterFile& context,
                                                    SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_getcwd))
    return std::nullopt;

  const address_t buffer_address = detail::syscall_arg(context, 0);
  const word_t size = detail::syscall_arg(context, 1);
  if (buffer_address == 0)
    return detail::return_error(context, detail::linux_efault);
  if (size == 0)
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, size))
    return detail::return_error(context, detail::linux_efault);

  try {
    std::string cwd = std::filesystem::current_path().string();
    if (cwd.size() + 1 > size)
      return detail::return_error(context, ERANGE);

    std::vector<std::byte> bytes(cwd.size() + 1);
    std::memcpy(bytes.data(), cwd.data(), cwd.size());
    auto written = machine.write_memory(buffer_address, bytes);
    if (!written)
      return detail::return_error(context, detail::memory_error_to_linux(written.error()));
    return detail::return_value(context, static_cast<std::int64_t>(cwd.size() + 1));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysTime::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_gettimeofday) &&
      !detail::handles(context, kind, detail::syscall_nanosleep) &&
      !detail::handles(context, kind, detail::syscall_getrusage) &&
      !detail::handles(context, kind, detail::syscall_clock_gettime)) {
    return std::nullopt;
  }

  const word_t number = detail::syscall_number(context);
  if (number == detail::syscall_nanosleep) {
    const address_t request_address = detail::syscall_arg(context, 0);
    if (request_address == 0)
      return detail::return_error(context, detail::linux_efault);

    std::array<std::byte, 16> bytes{};
    if (auto read_error = detail::read_guest_memory(machine, request_address, bytes))
      return detail::return_error(context, *read_error);

    timespec request{};
    request.tv_sec = static_cast<time_t>(static_cast<std::int64_t>(detail::read_le(bytes, 0, 8)));
    request.tv_nsec = static_cast<long>(detail::read_le(bytes, 8, 8));
    if (request.tv_sec < 0 || request.tv_nsec < 0 || request.tv_nsec >= 1000000000L)
      return detail::return_error(context, detail::linux_einval);

    if (nanosleep(&request, nullptr) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));
    return detail::return_value(context, 0);
  }

  if (number == detail::syscall_getrusage) {
    const address_t usage_address = detail::syscall_arg(context, 1);
    if (usage_address == 0)
      return detail::return_value(context, 0);

    std::array<std::byte, 144> bytes{};
    (void)machine.write_memory(usage_address, bytes);
    return detail::return_value(context, 0);
  }

  if (number == detail::syscall_gettimeofday) {
    timeval tv{};
    if (gettimeofday(&tv, nullptr) < 0)
      return detail::return_error(context, detail::host_errno_to_linux(errno));

    const address_t timeval_address = detail::syscall_arg(context, 0);
    if (timeval_address != 0) {
      std::array<std::byte, 16> bytes{};
      detail::write_le(bytes, 0, static_cast<std::uint64_t>(tv.tv_sec), 8);
      detail::write_le(bytes, 8, static_cast<std::uint64_t>(tv.tv_usec), 8);
      auto written = machine.write_memory(timeval_address, bytes);
      if (!written)
        return detail::return_error(context, detail::memory_error_to_linux(written.error()));
    }

    const address_t timezone_address = detail::syscall_arg(context, 1);
    if (timezone_address != 0) {
      std::array<std::byte, 8> bytes{};
      auto written = machine.write_memory(timezone_address, bytes);
      if (!written)
        return detail::return_error(context, detail::memory_error_to_linux(written.error()));
    }
    return detail::return_value(context, 0);
  }

  const address_t timespec_address = detail::syscall_arg(context, 1);
  if (timespec_address == 0)
    return detail::return_error(context, detail::linux_efault);

  timespec ts{};
  if (clock_gettime(CLOCK_REALTIME, &ts) < 0)
    return detail::return_error(context, detail::host_errno_to_linux(errno));

  std::array<std::byte, 16> bytes{};
  detail::write_le(bytes, 0, static_cast<std::uint64_t>(ts.tv_sec), 8);
  detail::write_le(bytes, 8, static_cast<std::uint64_t>(ts.tv_nsec), 8);
  auto written = machine.write_memory(timespec_address, bytes);
  if (!written)
    return detail::return_error(context, detail::memory_error_to_linux(written.error()));
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysSetTidAddress::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_set_tid_address))
    return std::nullopt;

  return detail::return_value(context, 1);
}

std::optional<SyscallResult> SysSetRobustList::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_set_robust_list))
    return std::nullopt;

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysRseq::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_rseq))
    return std::nullopt;

  return detail::return_error(context, detail::linux_enosys);
}

std::optional<SyscallResult> SysPrlimit64::try_syscall(Machine& machine, RegisterFile& context,
                                                       SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_prlimit64))
    return std::nullopt;

  const word_t pid = detail::syscall_arg(context, 0);
  if (pid != 0 && pid != detail::synthetic_pid)
    return detail::return_error(context, detail::linux_esrch);

  const word_t resource = detail::syscall_arg(context, 1);
  if (resource >= limits.size())
    return detail::return_error(context, detail::linux_einval);

  const address_t new_limit_address = detail::syscall_arg(context, 2);
  const address_t old_limit_address = detail::syscall_arg(context, 3);
  const Limit old_limit = limits[resource];
  Limit new_limit = old_limit;

  if (new_limit_address != 0) {
    if (detail::range_overflows(new_limit_address, 16))
      return detail::return_error(context, detail::linux_efault);

    std::array<std::byte, 16> bytes{};
    if (auto read_error = detail::read_guest_memory(machine, new_limit_address, bytes))
      return detail::return_error(context, *read_error);

    new_limit = detail::read_rlimit64(bytes);
    if (new_limit.current > new_limit.maximum)
      return detail::return_error(context, detail::linux_einval);
  }

  if (old_limit_address != 0) {
    if (detail::range_overflows(old_limit_address, 16))
      return detail::return_error(context, detail::linux_efault);
    if (auto error = detail::write_rlimit64(machine, old_limit_address, old_limit))
      return detail::return_error(context, *error);
  }

  if (new_limit_address != 0)
    limits[resource] = new_limit;

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysExit::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_exit))
    return std::nullopt;

  const int status = static_cast<int>(detail::syscall_arg(context, 0) & 0xff);
  if (exit_status)
    *exit_status = status;
  return SyscallResult{.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};
}

std::optional<SyscallResult> SysExitGroup::try_syscall(Machine&, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_exit_group))
    return std::nullopt;

  const int status = static_cast<int>(detail::syscall_arg(context, 0) & 0xff);
  if (exit_status)
    *exit_status = status;
  return SyscallResult{.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};
}

std::optional<SyscallResult> SysMmap::try_syscall(Machine& machine, RegisterFile& context, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_mmap))
    return std::nullopt;

  const address_t requested_address = detail::syscall_arg(context, 0);
  const word_t length = detail::syscall_arg(context, 1);
  const word_t prot = detail::syscall_arg(context, 2);
  const word_t flags = detail::syscall_arg(context, 3);
  const word_t offset = detail::syscall_arg(context, 5);

  if (length == 0 || !detail::page_aligned(offset))
    return detail::return_error(context, detail::linux_einval);

  const bool private_mapping = (flags & detail::linux_map_private) != 0;
  const bool shared_mapping = (flags & detail::linux_map_shared) != 0;
  if (private_mapping == shared_mapping)
    return detail::return_error(context, detail::linux_einval);

  constexpr word_t supported_flags = detail::linux_map_shared | detail::linux_map_private | detail::linux_map_fixed |
                                     detail::linux_map_anonymous | detail::linux_map_denied_write |
                                     detail::linux_map_executable | detail::linux_map_noreserve |
                                     detail::linux_map_populate | detail::linux_map_stack;
  if ((flags & ~supported_flags) != 0)
    return detail::return_error(context, detail::linux_einval);

  auto protection = detail::protection_from_linux(prot);
  if (!protection)
    return detail::return_error(context, detail::linux_einval);

  auto aligned_length = detail::page_align_up(length);
  if (!aligned_length)
    return detail::return_error(context, detail::linux_enomem);

  address_t mapped_address = 0;
  if ((flags & detail::linux_map_fixed) != 0) {
    if (!detail::page_aligned(requested_address))
      return detail::return_error(context, detail::linux_einval);
    mapped_address = requested_address;
  } else {
    auto next = detail::page_align_up(next_mapping_address);
    if (!next || *next > std::numeric_limits<address_t>::max() - *aligned_length)
      return detail::return_error(context, detail::linux_enomem);
    mapped_address = *next;
    next_mapping_address = mapped_address + *aligned_length;
  }

  if (mapped_address > std::numeric_limits<address_t>::max() - *aligned_length)
    return detail::return_error(context, detail::linux_einval);

  auto mapped = machine.map(mapped_address, *aligned_length, *protection);
  if (!mapped)
    return detail::return_error(context, detail::memory_error_to_linux(mapped.error()));

  if ((flags & detail::linux_map_anonymous) == 0) {
    const auto fd = detail::checked_fd(detail::syscall_arg(context, 4));
    if (!fd) {
      machine.unmap(mapped_address, *aligned_length);
      return detail::return_error(context, detail::linux_ebadf);
    }

    try {
      std::array<std::byte, Machine::kPageSize> buffer{};
      word_t total = 0;
      while (total < length) {
        const auto chunk_size =
            static_cast<std::size_t>(std::min<word_t>(length - total, static_cast<word_t>(buffer.size())));
        const ssize_t bytes_read = pread(*fd, buffer.data(), chunk_size, static_cast<off_t>(offset + total));
        if (bytes_read < 0) {
          const int saved_errno = errno;
          machine.unmap(mapped_address, *aligned_length);
          return detail::return_error(context, detail::host_errno_to_linux(saved_errno));
        }
        if (bytes_read == 0)
          break;

        auto written = machine.write_memory(
            mapped_address + total, std::span<const std::byte>(buffer.data(), static_cast<std::size_t>(bytes_read)));
        if (!written) {
          machine.unmap(mapped_address, *aligned_length);
          return detail::return_error(context, detail::memory_error_to_linux(written.error()));
        }
        total += static_cast<word_t>(bytes_read);
      }
    } catch (...) {
      machine.unmap(mapped_address, *aligned_length);
      return detail::return_error(context, detail::linux_eio);
    }
  }

  return detail::return_value(context, static_cast<std::int64_t>(mapped_address));
}

std::optional<SyscallResult> SysMunmap::try_syscall(Machine& machine, RegisterFile& context,
                                                    SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_munmap))
    return std::nullopt;

  const address_t address = detail::syscall_arg(context, 0);
  const word_t length = detail::syscall_arg(context, 1);
  if (length == 0 || !detail::page_aligned(address))
    return detail::return_error(context, detail::linux_einval);

  auto aligned_length = detail::page_align_up(length);
  if (!aligned_length)
    return detail::return_error(context, detail::linux_einval);

  machine.unmap(address, *aligned_length);
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysMremap::try_syscall(Machine& machine, RegisterFile& context,
                                                    SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_mremap))
    return std::nullopt;

  const address_t old_address = detail::syscall_arg(context, 0);
  const word_t old_size = detail::syscall_arg(context, 1);
  const word_t new_size = detail::syscall_arg(context, 2);
  const word_t flags = detail::syscall_arg(context, 3);
  const address_t fixed_address = detail::syscall_arg(context, 4);

  constexpr word_t supported_flags =
      detail::linux_mremap_maymove | detail::linux_mremap_fixed | detail::linux_mremap_dontunmap;
  if ((flags & ~supported_flags) != 0)
    return detail::return_error(context, detail::linux_einval);
  if ((flags & detail::linux_mremap_fixed) != 0 && (flags & detail::linux_mremap_maymove) == 0)
    return detail::return_error(context, detail::linux_einval);
  if ((flags & detail::linux_mremap_dontunmap) != 0)
    return detail::return_error(context, detail::linux_einval);
  if (old_size == 0 || new_size == 0 || !detail::page_aligned(old_address))
    return detail::return_error(context, detail::linux_einval);

  auto old_length = detail::page_align_up(old_size);
  auto new_length = detail::page_align_up(new_size);
  if (!old_length || !new_length)
    return detail::return_error(context, detail::linux_enomem);
  if (detail::range_overflows(old_address, *old_length))
    return detail::return_error(context, detail::linux_efault);

  try {
    const word_t copy_length = std::min(*old_length, *new_length);
    std::vector<std::byte> copied(static_cast<std::size_t>(copy_length));
    if (copy_length != 0) {
      auto read = machine.read_memory_into(old_address, copied);
      if (!read)
        return detail::return_error(context, detail::memory_error_to_linux(read.error()));
    }

    if (*new_length <= *old_length) {
      machine.unmap(old_address + *new_length, *old_length - *new_length);
      return detail::return_value(context, static_cast<std::int64_t>(old_address));
    }

    const address_t old_end = old_address + *old_length;
    const word_t growth = *new_length - *old_length;
    if (old_end <= std::numeric_limits<address_t>::max() - growth &&
        detail::range_is_unmapped(machine, old_end, growth)) {
      auto mapped = machine.map(old_end, growth, Protection::read | Protection::write);
      if (mapped)
        return detail::return_value(context, static_cast<std::int64_t>(old_address));
    }

    if ((flags & detail::linux_mremap_maymove) == 0)
      return detail::return_error(context, detail::linux_enomem);

    address_t new_address = 0;
    if ((flags & detail::linux_mremap_fixed) != 0) {
      if (!detail::page_aligned(fixed_address))
        return detail::return_error(context, detail::linux_einval);
      if (!detail::range_is_unmapped(machine, fixed_address, *new_length))
        machine.unmap(fixed_address, *new_length);
      new_address = fixed_address;
    } else {
      auto allocated = detail::find_unmapped_range(machine, next_mapping_address, *new_length);
      if (!allocated)
        return detail::return_error(context, detail::linux_enomem);
      new_address = *allocated;
    }

    auto mapped = machine.map(new_address, *new_length, Protection::read | Protection::write);
    if (!mapped)
      return detail::return_error(context, detail::memory_error_to_linux(mapped.error()));

    if (!copied.empty()) {
      auto written = machine.write_memory(new_address, copied);
      if (!written) {
        machine.unmap(new_address, *new_length);
        return detail::return_error(context, detail::memory_error_to_linux(written.error()));
      }
    }

    machine.unmap(old_address, *old_length);
    return detail::return_value(context, static_cast<std::int64_t>(new_address));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

} // namespace x86sim::linux_syscalls
