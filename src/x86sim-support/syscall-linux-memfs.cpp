// Memfs syscall personality: Linux file syscalls against the portable
// in-memory filesystem. Like syscall-linux.cpp this translation unit is fully
// portable — no host header may be included; every operation works on the
// shared MemFsState and the guest AddressSpace only.
#include "x86sim-support/syscall-linux-memfs.hpp"

#include "syscall-linux-detail.hpp"
#include "x86sim/addrspace.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <mutex>
#include <span>
#include <utility>
#include <vector>

namespace x86sim::linux_syscalls {

namespace {

using memfs::Node;

// Linux caps a single read/write transfer at MAX_RW_COUNT.
inline constexpr word_t max_rw_count = 0x7ffff000;

// getdents64 d_type values.
inline constexpr std::byte dt_dir{4};
inline constexpr std::byte dt_reg{8};
inline constexpr std::byte dt_lnk{10};

// utimensat tv_nsec sentinels.
inline constexpr std::uint64_t utime_now = 0x3fffffff;
inline constexpr std::uint64_t utime_omit = 0x3ffffffe;

[[nodiscard]] std::byte dirent_type(const Node& node) noexcept {
  if (node.is_dir())
    return dt_dir;
  if (node.is_symlink())
    return dt_lnk;
  return dt_reg;
}

// --- file-backed mmap loader shim ------------------------------------------
// The portable SysMmap delegates file-backed mappings to a single global
// function pointer. The shim installed here consults a pid -> state registry
// so it can serve memfs fds, and delegates anything else to whichever loader
// was registered before it (the host-POSIX one, or none). Registration is
// lazy — nothing changes for builds/machines that never enable memfs.
std::unordered_map<ProcessId, std::weak_ptr<MemFsState>>& mmap_registry() {
  static std::unordered_map<ProcessId, std::weak_ptr<MemFsState>> registry;
  return registry;
}

detail::MmapFileLoader g_previous_mmap_loader = nullptr;

std::optional<int> memfs_mmap_loader(ProcessId pid, AddressSpace& space, address_t dest, word_t length, word_t fd_arg,
                                     word_t offset) noexcept {
  std::shared_ptr<MemFsState> state;
  if (auto it = mmap_registry().find(pid); it != mmap_registry().end())
    state = it->second.lock();

  if (state) {
    if (auto proc_it = state->procs.find(pid); proc_it != state->procs.end()) {
      if (auto fd = detail::checked_fd(fd_arg)) {
        if (auto fd_it = proc_it->second.fds.find(*fd); fd_it != proc_it->second.fds.end()) {
          const auto& file = *fd_it->second.file;
          if (!file.node->is_file())
            return detail::linux_eacces;
          // Copy in bounded chunks; anything past EOF stays zero, matching the
          // zero-fill of a fresh mapping.
          std::array<std::byte, detail::io_chunk_size> chunk;
          word_t copied = 0;
          while (copied < length) {
            const std::size_t want =
                static_cast<std::size_t>(std::min<word_t>(length - copied, static_cast<word_t>(chunk.size())));
            const std::size_t got = state->fs->read(*file.node, offset + copied, std::span(chunk.data(), want));
            if (got == 0)
              break;
            auto written = space.write(dest + copied, std::span(chunk.data(), got));
            if (!written)
              return detail::memory_error_to_linux(written.error());
            copied += static_cast<word_t>(got);
            if (got < want)
              break;
          }
          return std::nullopt;
        }
      }
    }
  }

  if (g_previous_mmap_loader)
    return g_previous_mmap_loader(pid, space, dest, length, fd_arg, offset);
  return detail::linux_ebadf;
}

void install_mmap_shim() {
  static std::once_flag installed;
  std::call_once(installed, [] {
    g_previous_mmap_loader = detail::mmap_file_loader();
    detail::set_mmap_file_loader(&memfs_mmap_loader);
  });
}

} // namespace

// --- MemFsState -------------------------------------------------------------

MemFsState::MemFsState(std::shared_ptr<memfs::Filesystem> fs)
    : fs(std::move(fs)), mount(make_mount_namespace(this->fs)) {}

MemFsState::~MemFsState() {
  // Scrub the mmap registry so a later state cannot alias our pids.
  auto& registry = mmap_registry();
  for (const auto& [pid, proc] : procs)
    if (auto it = registry.find(pid); it != registry.end() && it->second.expired())
      registry.erase(it);
}

MemFsState::ProcState& MemFsState::proc(ProcessId pid) {
  auto [it, inserted] = procs.try_emplace(pid);
  if (inserted) {
    it->second.ns.mount = mount;
    mmap_registry()[pid] = weak_from_this();
  }
  return it->second;
}

void MemFsState::clone_proc(ProcessId parent, ProcessId child) {
  // Copying the fd table shares the OpenFile descriptions (shared offsets,
  // like fork()); namespaces are shared by shared_ptr copy.
  ProcState snapshot = proc(parent);
  procs.insert_or_assign(child, std::move(snapshot));
  mmap_registry()[child] = weak_from_this();
}

void MemFsState::erase_proc(ProcessId pid) {
  procs.erase(pid);
  if (auto it = mmap_registry().find(pid); it != mmap_registry().end() && it->second.expired())
    mmap_registry().erase(it);
}

int MemFsState::allocate_fd(const ProcState& proc, int minimum) const {
  int fd = std::max(minimum, 0);
  while (proc.fds.contains(fd) || (reserved_fds.contains(fd) && !proc.closed_reserved.contains(fd)))
    ++fd;
  return fd;
}

void MemFsState::install_at(ProcessId pid, int fd, OpenFile file, bool cloexec) {
  ProcState& state = proc(pid);
  state.fds.insert_or_assign(fd, FdEntry{.file = std::make_shared<OpenFile>(std::move(file)), .cloexec = cloexec});
  state.closed_reserved.erase(fd);
}

std::shared_ptr<MemFsState> make_memfs_state(std::shared_ptr<memfs::Filesystem> fs) {
  if (!fs)
    fs = std::make_shared<memfs::Filesystem>();
  install_mmap_shim();
  return std::make_shared<MemFsState>(std::move(fs));
}

// --- syscall handling --------------------------------------------------------

namespace {

namespace d = detail;

using ProcState = MemFsState::ProcState;
using OpenFile = MemFsState::OpenFile;
using FdEntry = MemFsState::FdEntry;

// fd classification: memfs-owned fds are handled here, tombstoned fds (a
// reserved fd the guest closed) fail with EBADF, anything else falls through
// to the rest of the chain.
enum class FdClass { ours, tombstone, foreign };

struct FdRef {
  FdClass kind = FdClass::foreign;
  FdEntry* entry = nullptr;
};

[[nodiscard]] FdRef classify_fd(ProcState& proc, word_t raw_fd) {
  const auto fd = d::checked_fd(raw_fd);
  if (!fd)
    return {FdClass::foreign, nullptr};
  if (auto it = proc.fds.find(*fd); it != proc.fds.end())
    return {FdClass::ours, &it->second};
  if (proc.closed_reserved.contains(*fd))
    return {FdClass::tombstone, nullptr};
  return {FdClass::foreign, nullptr};
}

// Base directory (a canonical path) for resolving `path` against `raw_dirfd`.
[[nodiscard]] std::expected<std::string, int> base_for(ProcState& proc, word_t raw_dirfd, std::string_view path) {
  if (!path.empty() && path.front() == '/')
    return std::string("/"); // absolute: base is ignored
  if (d::is_linux_at_fdcwd(raw_dirfd))
    return proc.cwd;
  const FdRef ref = classify_fd(proc, raw_dirfd);
  if (ref.kind != FdClass::ours)
    return std::unexpected(d::linux_ebadf);
  if (!ref.entry->file->node->is_dir())
    return std::unexpected(d::linux_enotdir);
  return ref.entry->file->path;
}

[[nodiscard]] std::expected<std::shared_ptr<Node>, int> resolve_at(ProcState& proc, word_t raw_dirfd,
                                                                   const std::string& path, bool follow) {
  auto base = base_for(proc, raw_dirfd, path);
  if (!base)
    return std::unexpected(base.error());
  const MountNamespace& mount = *proc.ns.mount;
  return mount.fs->resolve(path, mount.root, *base, follow);
}

[[nodiscard]] std::expected<memfs::Filesystem::ParentRef, int> parent_at(ProcState& proc, word_t raw_dirfd,
                                                                         const std::string& path) {
  auto base = base_for(proc, raw_dirfd, path);
  if (!base)
    return std::unexpected(base.error());
  const MountNamespace& mount = *proc.ns.mount;
  return mount.fs->resolve_parent(path, mount.root, *base);
}

[[nodiscard]] std::optional<int> write_memfs_stat(AddressSpace& space, address_t address, const memfs::StatInfo& info) {
  std::array<std::byte, 144> bytes{};
  d::write_le(bytes, 0, 1, 8); // st_dev: single synthetic device
  d::write_le(bytes, 8, info.ino, 8);
  d::write_le(bytes, 16, info.nlink, 8);
  d::write_le(bytes, 24, info.mode, 4);
  d::write_le(bytes, 28, d::synthetic_uid, 4);
  d::write_le(bytes, 32, d::synthetic_gid, 4);
  d::write_le(bytes, 40, 0, 8); // st_rdev
  d::write_le(bytes, 48, info.size, 8);
  d::write_le(bytes, 56, 4096, 8); // st_blksize
  d::write_le(bytes, 64, (info.size + 511) / 512, 8);
  d::write_le(bytes, 72, info.atime, 8);
  d::write_le(bytes, 80, 0, 8);
  d::write_le(bytes, 88, info.mtime, 8);
  d::write_le(bytes, 96, 0, 8);
  d::write_le(bytes, 104, info.ctime, 8);
  d::write_le(bytes, 112, 0, 8);

  auto written = space.write(address, bytes);
  if (!written)
    return d::memory_error_to_linux(written.error());
  return std::nullopt;
}

// --- open -------------------------------------------------------------------

[[nodiscard]] SyscallResult do_open(MemFsState& state, ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_openat = d::syscall_number(context) == d::syscall_openat;
  const word_t raw_dirfd = is_openat ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  const address_t path_address = d::syscall_arg(context, is_openat ? 1 : 0);
  const word_t flags = d::syscall_arg(context, is_openat ? 2 : 1);
  const word_t create_mode = d::syscall_arg(context, is_openat ? 3 : 2);

  constexpr word_t supported_flags = d::linux_o_accmode | d::linux_o_creat | d::linux_o_excl | d::linux_o_noctty |
                                     d::linux_o_trunc | d::linux_o_append | d::linux_o_nonblock | d::linux_o_dsync |
                                     d::linux_o_largefile | d::linux_o_directory | d::linux_o_nofollow |
                                     d::linux_o_cloexec | d::linux_o_sync;
  if ((flags & ~supported_flags) != 0)
    return d::return_error(context, d::linux_einval);
  const word_t accmode = flags & d::linux_o_accmode;
  if (accmode == d::linux_o_accmode)
    return d::return_error(context, d::linux_einval);

  auto path = d::read_c_string(space, path_address);
  if (!path.ok)
    return d::return_error(context, path.error);

  auto base = base_for(proc, raw_dirfd, path.value);
  if (!base)
    return d::return_error(context, base.error());
  const MountNamespace& mount = *proc.ns.mount;

  const bool follow = (flags & d::linux_o_nofollow) == 0;
  const bool wants_directory = (flags & d::linux_o_directory) != 0;
  const bool trailing_slash = !path.value.empty() && path.value.back() == '/';

  std::shared_ptr<Node> node;
  auto resolved = mount.fs->resolve(path.value, mount.root, *base, follow);
  if (resolved) {
    if ((flags & d::linux_o_creat) != 0 && (flags & d::linux_o_excl) != 0)
      return d::return_error(context, d::linux_eexist);
    node = *resolved;
  } else if (resolved.error() == d::linux_enoent && (flags & d::linux_o_creat) != 0) {
    if (trailing_slash)
      return d::return_error(context, d::linux_eisdir);
    auto parent = mount.fs->resolve_parent(path.value, mount.root, *base);
    if (!parent)
      return d::return_error(context, parent.error());
    const auto permissions = static_cast<std::uint32_t>(create_mode & 0777 & ~proc.umask);
    auto created = mount.fs->create_file(*parent->dir, parent->leaf, permissions);
    if (!created)
      return d::return_error(context, created.error());
    node = *created;
  } else {
    return d::return_error(context, resolved.error());
  }

  if (!follow && node->is_symlink())
    return d::return_error(context, d::linux_eloop);
  if (node->is_dir()) {
    if (accmode != d::linux_o_rdonly)
      return d::return_error(context, d::linux_eisdir);
  } else {
    if (wants_directory || trailing_slash)
      return d::return_error(context, d::linux_enotdir);
  }

  if ((flags & d::linux_o_trunc) != 0 && node->is_file() && accmode != d::linux_o_rdonly) {
    if (auto truncated = mount.fs->truncate(*node, 0); !truncated)
      return d::return_error(context, truncated.error());
  }

  OpenFile file{.node = node, .path = {}, .offset = 0, .status_flags = flags & ~(d::linux_o_creat | d::linux_o_excl |
                                                                                 d::linux_o_trunc | d::linux_o_cloexec)};
  if (node->is_dir()) {
    auto canonical = mount.fs->canonicalize_directory(path.value, mount.root, *base);
    if (!canonical)
      return d::return_error(context, canonical.error());
    file.path = std::move(*canonical);
  }

  const int fd = state.allocate_fd(proc);
  proc.fds.insert_or_assign(fd, FdEntry{.file = std::make_shared<OpenFile>(std::move(file)),
                                        .cloexec = (flags & d::linux_o_cloexec) != 0});
  proc.closed_reserved.erase(fd);
  return d::return_value(context, fd);
}

// --- read / write -----------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_read(MemFsState& state, ProcState& proc, CpuState& context,
                                                   AddressSpace& space, bool positional) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  OpenFile& file = *ref.entry->file;
  if (file.node->is_dir())
    return d::return_error(context, d::linux_eisdir);
  if ((file.status_flags & d::linux_o_accmode) == d::linux_o_wronly)
    return d::return_error(context, d::linux_ebadf);

  const address_t buffer = d::syscall_arg(context, 1);
  const word_t count = std::min(d::syscall_arg(context, 2), max_rw_count);
  const std::uint64_t offset = positional ? d::syscall_arg(context, 3) : file.offset;
  if (positional && static_cast<std::int64_t>(offset) < 0)
    return d::return_error(context, d::linux_einval);
  if (count == 0)
    return d::return_value(context, 0);
  if (buffer == 0 || d::range_overflows(buffer, count))
    return d::return_error(context, d::linux_efault);

  std::vector<std::byte> data(static_cast<std::size_t>(count));
  const std::size_t got = state.fs->read(*file.node, offset, data);
  if (got > 0) {
    auto written = space.write(buffer, std::span(data.data(), got));
    if (!written)
      return d::return_error(context, d::memory_error_to_linux(written.error()));
  }
  if (!positional)
    file.offset = offset + got;
  return d::return_value(context, static_cast<std::int64_t>(got));
}

[[nodiscard]] std::optional<SyscallResult> do_write(MemFsState& state, ProcState& proc, CpuState& context,
                                                    AddressSpace& space, bool positional) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  OpenFile& file = *ref.entry->file;
  if ((file.status_flags & d::linux_o_accmode) == d::linux_o_rdonly)
    return d::return_error(context, d::linux_ebadf);

  const address_t buffer = d::syscall_arg(context, 1);
  const word_t count = std::min(d::syscall_arg(context, 2), max_rw_count);
  std::uint64_t offset = positional ? d::syscall_arg(context, 3) : file.offset;
  if (positional && static_cast<std::int64_t>(offset) < 0)
    return d::return_error(context, d::linux_einval);
  // O_APPEND writes atomically at end-of-file (pwrite on Linux ignores the
  // offset for O_APPEND files as well).
  if ((file.status_flags & d::linux_o_append) != 0)
    offset = file.node->size();
  if (count == 0)
    return d::return_value(context, 0);
  if (buffer == 0 || d::range_overflows(buffer, count))
    return d::return_error(context, d::linux_efault);

  std::vector<std::byte> data(static_cast<std::size_t>(count));
  auto read_guest = space.read(buffer, std::span(data));
  if (!read_guest)
    return d::return_error(context, d::memory_error_to_linux(read_guest.error()));

  auto written = state.fs->write(*file.node, offset, data);
  if (!written)
    return d::return_error(context, written.error());
  if (!positional)
    file.offset = offset + *written;
  return d::return_value(context, static_cast<std::int64_t>(*written));
}

// --- close / lseek ----------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_close(MemFsState& state, ProcState& proc, CpuState& context) {
  const auto fd = d::checked_fd(d::syscall_arg(context, 0));
  if (!fd)
    return d::return_error(context, d::linux_ebadf);

  if (auto it = proc.fds.find(*fd); it != proc.fds.end()) {
    proc.fds.erase(it);
    // Closing a memfs fd that shadowed an embedder fd must not resurrect the
    // embedder's stream: the number stays dead until reused.
    if (state.reserved_fds.contains(*fd))
      proc.closed_reserved.insert(*fd);
    return d::return_value(context, 0);
  }
  if (state.reserved_fds.contains(*fd) && !proc.closed_reserved.contains(*fd)) {
    // Closing an embedder-owned fd: the stream stays alive on the host side,
    // but the guest fd number is dead from now on.
    proc.closed_reserved.insert(*fd);
    return d::return_value(context, 0);
  }
  return d::return_error(context, d::linux_ebadf);
}

[[nodiscard]] std::optional<SyscallResult> do_lseek(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  OpenFile& file = *ref.entry->file;
  const auto offset = static_cast<std::int64_t>(d::syscall_arg(context, 1));
  const word_t whence = d::syscall_arg(context, 2);

  if (file.node->is_dir()) {
    // Directory offsets are opaque cookies; only rewinddir's lseek(fd, 0,
    // SEEK_SET) is supported.
    if (whence == 0 && offset == 0) {
      file.dir_cursor.clear();
      file.dir_phase = 0;
      file.dir_cookie = 0;
      file.offset = 0;
      return d::return_value(context, 0);
    }
    if (whence == 1 && offset == 0)
      return d::return_value(context, static_cast<std::int64_t>(file.dir_cookie));
    return d::return_error(context, d::linux_einval);
  }

  std::int64_t base = 0;
  switch (whence) {
  case 0: // SEEK_SET
    base = 0;
    break;
  case 1: // SEEK_CUR
    base = static_cast<std::int64_t>(file.offset);
    break;
  case 2: // SEEK_END
    base = static_cast<std::int64_t>(file.node->size());
    break;
  default:
    return d::return_error(context, d::linux_einval);
  }
  const std::int64_t target = base + offset;
  if (target < 0)
    return d::return_error(context, d::linux_einval);
  file.offset = static_cast<std::uint64_t>(target);
  return d::return_value(context, target);
}

// --- stat family --------------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_fstat(MemFsState& state, ProcState& proc, CpuState& context,
                                                    AddressSpace& space) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  if (auto error = write_memfs_stat(space, d::syscall_arg(context, 1), state.fs->stat(*ref.entry->file->node)))
    return d::return_error(context, *error);
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_stat(MemFsState& state, ProcState& proc, CpuState& context, AddressSpace& space) {
  const word_t number = d::syscall_number(context);
  const bool is_newfstatat = number == d::syscall_newfstatat;
  const word_t raw_dirfd = is_newfstatat ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  const address_t path_address = d::syscall_arg(context, is_newfstatat ? 1 : 0);
  const address_t stat_address = d::syscall_arg(context, is_newfstatat ? 2 : 1);
  const word_t flags = is_newfstatat ? d::syscall_arg(context, 3) : 0;

  constexpr word_t supported_flags = d::linux_at_symlink_nofollow | d::linux_at_empty_path;
  if ((flags & ~supported_flags) != 0)
    return d::return_error(context, d::linux_einval);
  const bool follow = number == d::syscall_lstat ? false : (flags & d::linux_at_symlink_nofollow) == 0;

  auto path = d::read_c_string(space, path_address);
  if (!path.ok)
    return d::return_error(context, path.error);

  std::shared_ptr<Node> node;
  if (is_newfstatat && (flags & d::linux_at_empty_path) != 0 && path.value.empty()) {
    if (d::is_linux_at_fdcwd(raw_dirfd)) {
      auto resolved = resolve_at(proc, d::linux_at_fdcwd, ".", true);
      if (!resolved)
        return d::return_error(context, resolved.error());
      node = *resolved;
    } else {
      const FdRef ref = classify_fd(proc, raw_dirfd);
      if (ref.kind != FdClass::ours)
        return d::return_error(context, d::linux_ebadf);
      node = ref.entry->file->node;
    }
  } else {
    auto resolved = resolve_at(proc, raw_dirfd, path.value, follow);
    if (!resolved)
      return d::return_error(context, resolved.error());
    node = *resolved;
  }

  if (auto error = write_memfs_stat(space, stat_address, state.fs->stat(*node)))
    return d::return_error(context, *error);
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_statx(MemFsState& state, ProcState& proc, CpuState& context, AddressSpace& space) {
  const word_t raw_dirfd = d::syscall_arg(context, 0);
  const address_t path_address = d::syscall_arg(context, 1);
  const word_t flags = d::syscall_arg(context, 2);
  const address_t statx_address = d::syscall_arg(context, 4);

  constexpr word_t supported_flags = d::linux_at_symlink_nofollow | d::linux_at_empty_path | 0x6000 /*sync flags*/;
  if ((flags & ~supported_flags) != 0)
    return d::return_error(context, d::linux_einval);
  const bool follow = (flags & d::linux_at_symlink_nofollow) == 0;

  auto path = d::read_c_string(space, path_address);
  if (!path.ok)
    return d::return_error(context, path.error);

  std::shared_ptr<Node> node;
  if ((flags & d::linux_at_empty_path) != 0 && path.value.empty()) {
    const FdRef ref = classify_fd(proc, raw_dirfd);
    if (ref.kind != FdClass::ours)
      return d::return_error(context, d::linux_ebadf);
    node = ref.entry->file->node;
  } else {
    auto resolved = resolve_at(proc, raw_dirfd, path.value, follow);
    if (!resolved)
      return d::return_error(context, resolved.error());
    node = *resolved;
  }

  const memfs::StatInfo info = state.fs->stat(*node);
  std::array<std::byte, 256> bytes{};
  d::write_le(bytes, 0, 0x7ff, 4); // stx_mask: STATX_BASIC_STATS
  d::write_le(bytes, 4, 4096, 4);  // stx_blksize
  d::write_le(bytes, 16, info.nlink, 4);
  d::write_le(bytes, 20, d::synthetic_uid, 4);
  d::write_le(bytes, 24, d::synthetic_gid, 4);
  d::write_le(bytes, 28, info.mode, 2);
  d::write_le(bytes, 32, info.ino, 8);
  d::write_le(bytes, 40, info.size, 8);
  d::write_le(bytes, 48, (info.size + 511) / 512, 8);
  d::write_le(bytes, 64, info.atime, 8);  // stx_atime.tv_sec
  d::write_le(bytes, 96, info.ctime, 8);  // stx_ctime.tv_sec
  d::write_le(bytes, 112, info.mtime, 8); // stx_mtime.tv_sec

  auto written = space.write(statx_address, bytes);
  if (!written)
    return d::return_error(context, d::memory_error_to_linux(written.error()));
  return d::return_value(context, 0);
}

// --- getdents64 ---------------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_getdents64(ProcState& proc, CpuState& context, AddressSpace& space) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  OpenFile& file = *ref.entry->file;
  if (!file.node->is_dir())
    return d::return_error(context, d::linux_enotdir);

  const address_t buffer_address = d::syscall_arg(context, 1);
  const word_t count = d::syscall_arg(context, 2);
  if (!d::fits_host_transfer(count))
    return d::return_error(context, d::linux_einval);
  if (buffer_address == 0 || d::range_overflows(buffer_address, count))
    return d::return_error(context, d::linux_efault);

  std::vector<std::byte> output;
  const auto& entries = file.node->dir().entries;

  auto emit = [&](const std::string& name, std::uint64_t ino, std::byte type) -> bool {
    const std::size_t record_length = (19 + name.size() + 1 + 7) & ~std::size_t{7};
    if (record_length > count - output.size())
      return false;
    const std::size_t offset = output.size();
    output.resize(offset + record_length, std::byte{0});
    auto write_field = [&](std::size_t field_offset, std::uint64_t value, std::size_t width) noexcept {
      for (std::size_t i = 0; i < width; ++i)
        output[offset + field_offset + i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);
    };
    write_field(0, ino, 8);
    write_field(8, ++file.dir_cookie, 8); // d_off: opaque, monotonically increasing
    write_field(16, record_length, 2);
    output[offset + 18] = type;
    std::memcpy(output.data() + offset + 19, name.data(), name.size());
    return true;
  };

  // "." and ".." are synthesized first (the parent's identity is not tracked,
  // so ".." reuses the directory's own inode — like a mount root).
  static const std::string dot = ".";
  static const std::string dotdot = "..";
  bool full = false;
  if (file.dir_phase == 0 && !full) {
    if (emit(dot, file.node->ino, dt_dir))
      file.dir_phase = 1;
    else
      full = true;
  }
  if (file.dir_phase == 1 && !full) {
    if (emit(dotdot, file.node->ino, dt_dir))
      file.dir_phase = 2;
    else
      full = true;
  }
  if (file.dir_phase == 2 && !full) {
    auto it = file.dir_cursor.empty() ? entries.begin() : entries.upper_bound(file.dir_cursor);
    for (; it != entries.end(); ++it) {
      if (!emit(it->first, it->second->ino, dirent_type(*it->second))) {
        full = true;
        break;
      }
      file.dir_cursor = it->first;
    }
  }

  if (output.empty()) {
    // Not even one record fits into a non-empty buffer request.
    if (full)
      return d::return_error(context, d::linux_einval);
    return d::return_value(context, 0);
  }

  auto written = space.write(buffer_address, output);
  if (!written)
    return d::return_error(context, d::memory_error_to_linux(written.error()));
  return d::return_value(context, static_cast<std::int64_t>(output.size()));
}

// --- directory / namespace mutations ------------------------------------------

[[nodiscard]] SyscallResult do_mkdir(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_mkdirat;
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const word_t mode = d::syscall_arg(context, is_at ? 2 : 1);

  auto parent = parent_at(proc, raw_dirfd, path.value);
  if (!parent)
    return d::return_error(context, parent.error());
  auto created = proc.ns.mount->fs->make_directory(*parent->dir, parent->leaf,
                                                   static_cast<std::uint32_t>(mode & 0777 & ~proc.umask));
  if (!created)
    return d::return_error(context, created.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_unlink(ProcState& proc, CpuState& context, AddressSpace& space) {
  const word_t number = d::syscall_number(context);
  const bool is_at = number == d::syscall_unlinkat;
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const word_t flags = is_at ? d::syscall_arg(context, 2) : 0;
  if ((flags & ~d::linux_at_removedir) != 0)
    return d::return_error(context, d::linux_einval);
  const bool remove_directory = number == d::syscall_rmdir || (flags & d::linux_at_removedir) != 0;

  auto parent = parent_at(proc, raw_dirfd, path.value);
  if (!parent)
    return d::return_error(context, parent.error());
  auto& fs = *proc.ns.mount->fs;
  auto removed = remove_directory ? fs.remove_directory(*parent->dir, parent->leaf)
                                  : fs.unlink(*parent->dir, parent->leaf);
  if (!removed)
    return d::return_error(context, removed.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_rename(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_renameat;
  const word_t old_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto old_path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!old_path.ok)
    return d::return_error(context, old_path.error);
  const word_t new_dirfd = is_at ? d::syscall_arg(context, 2) : d::linux_at_fdcwd;
  auto new_path = d::read_c_string(space, d::syscall_arg(context, is_at ? 3 : 1));
  if (!new_path.ok)
    return d::return_error(context, new_path.error);

  auto old_parent = parent_at(proc, old_dirfd, old_path.value);
  if (!old_parent)
    return d::return_error(context, old_parent.error());
  auto new_parent = parent_at(proc, new_dirfd, new_path.value);
  if (!new_parent)
    return d::return_error(context, new_parent.error());

  auto renamed =
      proc.ns.mount->fs->rename(*old_parent->dir, old_parent->leaf, *new_parent->dir, new_parent->leaf);
  if (!renamed)
    return d::return_error(context, renamed.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_symlink(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_symlinkat;
  auto target = d::read_c_string(space, d::syscall_arg(context, 0));
  if (!target.ok)
    return d::return_error(context, target.error);
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 1) : d::linux_at_fdcwd;
  auto link_path = d::read_c_string(space, d::syscall_arg(context, is_at ? 2 : 1));
  if (!link_path.ok)
    return d::return_error(context, link_path.error);

  auto parent = parent_at(proc, raw_dirfd, link_path.value);
  if (!parent)
    return d::return_error(context, parent.error());
  auto created = proc.ns.mount->fs->make_symlink(*parent->dir, parent->leaf, std::move(target.value));
  if (!created)
    return d::return_error(context, created.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_link(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_linkat;
  const word_t old_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto old_path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!old_path.ok)
    return d::return_error(context, old_path.error);
  const word_t new_dirfd = is_at ? d::syscall_arg(context, 2) : d::linux_at_fdcwd;
  auto new_path = d::read_c_string(space, d::syscall_arg(context, is_at ? 3 : 1));
  if (!new_path.ok)
    return d::return_error(context, new_path.error);
  const word_t flags = is_at ? d::syscall_arg(context, 4) : 0;
  constexpr word_t at_symlink_follow = 0x400;
  if ((flags & ~at_symlink_follow) != 0)
    return d::return_error(context, d::linux_einval);

  auto node = resolve_at(proc, old_dirfd, old_path.value, (flags & at_symlink_follow) != 0);
  if (!node)
    return d::return_error(context, node.error());
  auto parent = parent_at(proc, new_dirfd, new_path.value);
  if (!parent)
    return d::return_error(context, parent.error());
  auto linked = proc.ns.mount->fs->link(*parent->dir, parent->leaf, *node);
  if (!linked)
    return d::return_error(context, linked.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_readlink(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_readlinkat;
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const address_t buffer = d::syscall_arg(context, is_at ? 2 : 1);
  const word_t size = d::syscall_arg(context, is_at ? 3 : 2);
  if (size == 0 || static_cast<std::int64_t>(size) < 0)
    return d::return_error(context, d::linux_einval);

  auto node = resolve_at(proc, raw_dirfd, path.value, /*follow=*/false);
  if (!node)
    return d::return_error(context, node.error());
  if (!(*node)->is_symlink())
    return d::return_error(context, d::linux_einval);

  const std::string& target = (*node)->symlink().target;
  // readlink(2) truncates silently and does not NUL-terminate.
  const std::size_t n = std::min<std::size_t>(target.size(), static_cast<std::size_t>(size));
  if (n > 0) {
    auto written = space.write(buffer, std::as_bytes(std::span(target.data(), n)));
    if (!written)
      return d::return_error(context, d::memory_error_to_linux(written.error()));
  }
  return d::return_value(context, static_cast<std::int64_t>(n));
}

// --- cwd ------------------------------------------------------------------------

[[nodiscard]] SyscallResult do_chdir(ProcState& proc, CpuState& context, AddressSpace& space) {
  auto path = d::read_c_string(space, d::syscall_arg(context, 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const MountNamespace& mount = *proc.ns.mount;
  auto canonical = mount.fs->canonicalize_directory(path.value, mount.root, proc.cwd);
  if (!canonical)
    return d::return_error(context, canonical.error());
  proc.cwd = std::move(*canonical);
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<SyscallResult> do_fchdir(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);
  if (!ref.entry->file->node->is_dir())
    return d::return_error(context, d::linux_enotdir);
  proc.cwd = ref.entry->file->path;
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_getcwd(ProcState& proc, CpuState& context, AddressSpace& space) {
  const address_t buffer = d::syscall_arg(context, 0);
  const word_t size = d::syscall_arg(context, 1);
  const std::size_t needed = proc.cwd.size() + 1;
  if (size < needed)
    return d::return_error(context, size == 0 ? d::linux_einval : d::linux_erange);
  if (buffer == 0)
    return d::return_error(context, d::linux_efault);

  std::vector<std::byte> bytes(needed);
  std::memcpy(bytes.data(), proc.cwd.data(), proc.cwd.size());
  bytes[needed - 1] = std::byte{0};
  auto written = space.write(buffer, bytes);
  if (!written)
    return d::return_error(context, d::memory_error_to_linux(written.error()));
  return d::return_value(context, static_cast<std::int64_t>(needed));
}

// --- truncate / permissions / times ----------------------------------------------

[[nodiscard]] SyscallResult do_truncate(ProcState& proc, CpuState& context, AddressSpace& space) {
  auto path = d::read_c_string(space, d::syscall_arg(context, 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const auto length = static_cast<std::int64_t>(d::syscall_arg(context, 1));
  if (length < 0)
    return d::return_error(context, d::linux_einval);

  auto node = resolve_at(proc, d::linux_at_fdcwd, path.value, true);
  if (!node)
    return d::return_error(context, node.error());
  auto truncated = proc.ns.mount->fs->truncate(**node, static_cast<std::uint64_t>(length));
  if (!truncated)
    return d::return_error(context, truncated.error());
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<SyscallResult> do_ftruncate(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  OpenFile& file = *ref.entry->file;
  if ((file.status_flags & d::linux_o_accmode) == d::linux_o_rdonly)
    return d::return_error(context, d::linux_einval);
  const auto length = static_cast<std::int64_t>(d::syscall_arg(context, 1));
  if (length < 0)
    return d::return_error(context, d::linux_einval);
  auto truncated = proc.ns.mount->fs->truncate(*file.node, static_cast<std::uint64_t>(length));
  if (!truncated)
    return d::return_error(context, truncated.error());
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_chmod_path(ProcState& proc, CpuState& context, AddressSpace& space) {
  const bool is_at = d::syscall_number(context) == d::syscall_fchmodat;
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const word_t mode = d::syscall_arg(context, is_at ? 2 : 1);

  auto node = resolve_at(proc, raw_dirfd, path.value, true);
  if (!node)
    return d::return_error(context, node.error());
  proc.ns.mount->fs->set_permissions(**node, static_cast<std::uint32_t>(mode));
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<SyscallResult> do_fchmod(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);
  proc.ns.mount->fs->set_permissions(*ref.entry->file->node, static_cast<std::uint32_t>(d::syscall_arg(context, 1)));
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_umask(ProcState& proc, CpuState& context) {
  const word_t previous = proc.umask;
  proc.umask = d::syscall_arg(context, 0) & 0777;
  return d::return_value(context, static_cast<std::int64_t>(previous));
}

[[nodiscard]] SyscallResult do_utime(ProcState& proc, CpuState& context, AddressSpace& space) {
  auto path = d::read_c_string(space, d::syscall_arg(context, 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  auto node = resolve_at(proc, d::linux_at_fdcwd, path.value, true);
  if (!node)
    return d::return_error(context, node.error());

  auto& fs = *proc.ns.mount->fs;
  const address_t times_address = d::syscall_arg(context, 1);
  if (times_address == 0) {
    fs.touch(**node);
    return d::return_value(context, 0);
  }
  std::array<std::byte, 16> times{}; // struct utimbuf: actime, modtime
  if (auto error = d::read_guest_memory(space, times_address, times))
    return d::return_error(context, *error);
  fs.set_times(**node, d::read_le(times, 0, 8), d::read_le(times, 8, 8));
  return d::return_value(context, 0);
}

[[nodiscard]] SyscallResult do_utimensat(ProcState& proc, CpuState& context, AddressSpace& space) {
  const word_t raw_dirfd = d::syscall_arg(context, 0);
  const address_t path_address = d::syscall_arg(context, 1);
  const address_t times_address = d::syscall_arg(context, 2);
  const word_t flags = d::syscall_arg(context, 3);
  if ((flags & ~d::linux_at_symlink_nofollow) != 0)
    return d::return_error(context, d::linux_einval);

  std::shared_ptr<Node> node;
  if (path_address == 0) {
    // futimens(): a NULL path operates on the dirfd itself.
    const FdRef ref = classify_fd(proc, raw_dirfd);
    if (ref.kind != FdClass::ours)
      return d::return_error(context, d::linux_ebadf);
    node = ref.entry->file->node;
  } else {
    auto path = d::read_c_string(space, path_address);
    if (!path.ok)
      return d::return_error(context, path.error);
    auto resolved = resolve_at(proc, raw_dirfd, path.value, (flags & d::linux_at_symlink_nofollow) == 0);
    if (!resolved)
      return d::return_error(context, resolved.error());
    node = *resolved;
  }

  auto& fs = *proc.ns.mount->fs;
  if (times_address == 0) {
    fs.touch(*node);
    return d::return_value(context, 0);
  }
  std::array<std::byte, 32> times{}; // struct timespec[2]
  if (auto error = d::read_guest_memory(space, times_address, times))
    return d::return_error(context, *error);
  const std::uint64_t atime_nsec = d::read_le(times, 8, 8);
  const std::uint64_t mtime_nsec = d::read_le(times, 24, 8);
  const std::uint64_t now = fs.current_tick();
  const std::uint64_t atime = atime_nsec == utime_now    ? now
                              : atime_nsec == utime_omit ? node->atime
                                                         : d::read_le(times, 0, 8);
  const std::uint64_t mtime = mtime_nsec == utime_now    ? now
                              : mtime_nsec == utime_omit ? node->mtime
                                                         : d::read_le(times, 16, 8);
  fs.set_times(*node, atime, mtime);
  return d::return_value(context, 0);
}

// --- access ------------------------------------------------------------------------

[[nodiscard]] SyscallResult do_access(ProcState& proc, CpuState& context, AddressSpace& space) {
  const word_t number = d::syscall_number(context);
  const bool is_at = number == d::syscall_faccessat || number == d::syscall_faccessat2;
  const word_t raw_dirfd = is_at ? d::syscall_arg(context, 0) : d::linux_at_fdcwd;
  auto path = d::read_c_string(space, d::syscall_arg(context, is_at ? 1 : 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  const word_t flags = number == d::syscall_faccessat2 ? d::syscall_arg(context, 3) : 0;

  // Permissions are stored but not enforced: existence is the whole check.
  auto node = resolve_at(proc, raw_dirfd, path.value, (flags & d::linux_at_symlink_nofollow) == 0);
  if (!node)
    return d::return_error(context, node.error());
  return d::return_value(context, 0);
}

// --- dup / fcntl / ioctl --------------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_dup(MemFsState& state, ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  const int fd = state.allocate_fd(proc);
  proc.fds.insert_or_assign(fd, FdEntry{.file = ref.entry->file, .cloexec = false});
  proc.closed_reserved.erase(fd);
  return d::return_value(context, fd);
}

[[nodiscard]] std::optional<SyscallResult> do_dup2(MemFsState& state, ProcState& proc, CpuState& context) {
  const bool is_dup3 = d::syscall_number(context) == d::syscall_dup3;
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  const auto new_fd = d::checked_fd(d::syscall_arg(context, 1));
  if (!new_fd)
    return d::return_error(context, d::linux_ebadf);
  const word_t flags = is_dup3 ? d::syscall_arg(context, 2) : 0;
  if ((flags & ~d::linux_o_cloexec) != 0)
    return d::return_error(context, d::linux_einval);

  const int old_fd = static_cast<int>(d::syscall_arg(context, 0));
  if (old_fd == *new_fd) {
    if (is_dup3)
      return d::return_error(context, d::linux_einval);
    return d::return_value(context, *new_fd);
  }

  // Installing over any number is fine: a memfs entry shadows an embedder
  // stream at the same number because SysMemFs is consulted first.
  proc.fds.insert_or_assign(*new_fd, FdEntry{.file = ref.entry->file, .cloexec = (flags & d::linux_o_cloexec) != 0});
  proc.closed_reserved.erase(*new_fd);
  return d::return_value(context, *new_fd);
}

[[nodiscard]] std::optional<SyscallResult> do_fcntl(MemFsState& state, ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  const word_t command = d::syscall_arg(context, 1);
  const word_t argument = d::syscall_arg(context, 2);
  switch (command) {
  case d::linux_f_dupfd:
  case d::linux_f_dupfd_cloexec: {
    const auto minimum = d::checked_fd(argument);
    if (!minimum)
      return d::return_error(context, d::linux_einval);
    const int fd = state.allocate_fd(proc, *minimum);
    proc.fds.insert_or_assign(fd,
                              FdEntry{.file = ref.entry->file, .cloexec = command == d::linux_f_dupfd_cloexec});
    proc.closed_reserved.erase(fd);
    return d::return_value(context, fd);
  }
  case d::linux_f_getfd:
    return d::return_value(context, ref.entry->cloexec ? d::linux_fd_cloexec : 0);
  case d::linux_f_setfd:
    ref.entry->cloexec = (argument & d::linux_fd_cloexec) != 0;
    return d::return_value(context, 0);
  case d::linux_f_getfl:
    return d::return_value(context, static_cast<std::int64_t>(ref.entry->file->status_flags));
  case d::linux_f_setfl: {
    constexpr word_t mutable_flags = d::linux_o_append | d::linux_o_nonblock | d::linux_o_dsync | d::linux_o_sync;
    OpenFile& file = *ref.entry->file;
    file.status_flags = (file.status_flags & ~mutable_flags) | (argument & mutable_flags);
    return d::return_value(context, 0);
  }
  default:
    return d::return_error(context, d::linux_einval);
  }
}

// --- no-op durability / statfs ----------------------------------------------------------

[[nodiscard]] std::optional<SyscallResult> do_fsync_like(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);
  return d::return_value(context, 0); // memory is as durable as it gets
}

[[nodiscard]] std::optional<SyscallResult> do_fallocate(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);

  const word_t mode = d::syscall_arg(context, 1);
  const auto offset = static_cast<std::int64_t>(d::syscall_arg(context, 2));
  const auto length = static_cast<std::int64_t>(d::syscall_arg(context, 3));
  if (mode != 0)
    return d::return_error(context, d::linux_enotty); // only plain allocation is modeled
  if (offset < 0 || length <= 0)
    return d::return_error(context, d::linux_einval);

  OpenFile& file = *ref.entry->file;
  if (!file.node->is_file())
    return d::return_error(context, d::linux_ebadf);
  const std::uint64_t end = static_cast<std::uint64_t>(offset) + static_cast<std::uint64_t>(length);
  if (end > file.node->size()) {
    auto grown = proc.ns.mount->fs->truncate(*file.node, end);
    if (!grown)
      return d::return_error(context, grown.error());
  }
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<int> write_memfs_statfs(const MemFsState& state, AddressSpace& space, address_t address) {
  std::array<std::byte, 120> bytes{};
  const std::uint64_t block_size = 4096;
  const std::uint64_t blocks = state.fs->max_file_size() / block_size;
  d::write_le(bytes, 0, 0x01021994, 8); // f_type: TMPFS_MAGIC
  d::write_le(bytes, 8, block_size, 8); // f_bsize
  d::write_le(bytes, 16, blocks, 8);    // f_blocks
  d::write_le(bytes, 24, blocks, 8);    // f_bfree
  d::write_le(bytes, 32, blocks, 8);    // f_bavail
  d::write_le(bytes, 40, 0, 8);         // f_files
  d::write_le(bytes, 48, 0, 8);         // f_ffree
  d::write_le(bytes, 72, 255, 8);       // f_namelen
  d::write_le(bytes, 80, block_size, 8); // f_frsize

  auto written = space.write(address, bytes);
  if (!written)
    return d::memory_error_to_linux(written.error());
  return std::nullopt;
}

[[nodiscard]] SyscallResult do_statfs(MemFsState& state, ProcState& proc, CpuState& context, AddressSpace& space) {
  auto path = d::read_c_string(space, d::syscall_arg(context, 0));
  if (!path.ok)
    return d::return_error(context, path.error);
  auto node = resolve_at(proc, d::linux_at_fdcwd, path.value, true);
  if (!node)
    return d::return_error(context, node.error());
  if (auto error = write_memfs_statfs(state, space, d::syscall_arg(context, 1)))
    return d::return_error(context, *error);
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<SyscallResult> do_fstatfs(MemFsState& state, ProcState& proc, CpuState& context,
                                                      AddressSpace& space) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);
  if (auto error = write_memfs_statfs(state, space, d::syscall_arg(context, 1)))
    return d::return_error(context, *error);
  return d::return_value(context, 0);
}

[[nodiscard]] std::optional<SyscallResult> do_ioctl(ProcState& proc, CpuState& context) {
  const FdRef ref = classify_fd(proc, d::syscall_arg(context, 0));
  if (ref.kind == FdClass::foreign)
    return std::nullopt;
  if (ref.kind == FdClass::tombstone)
    return d::return_error(context, d::linux_ebadf);
  return d::return_error(context, d::linux_enotty); // memfs files are never terminals
}

} // namespace

std::optional<SyscallResult> SysMemFs::try_syscall(Machine&, ProcessId pid, CpuState& context, AddressSpace& space,
                                                   SyscallKind kind) noexcept {
  if (kind != SyscallKind::syscall64)
    return std::nullopt;
  if (!state)
    return std::nullopt;

  const word_t number = d::syscall_number(context);
  try {
    ProcState& proc = state->proc(pid);

    switch (number) {
    // fd-based: fall through (nullopt) for fds that are not ours.
    case d::syscall_read:
      return do_read(*state, proc, context, space, /*positional=*/false);
    case d::syscall_pread64:
      return do_read(*state, proc, context, space, /*positional=*/true);
    case d::syscall_write:
      return do_write(*state, proc, context, space, /*positional=*/false);
    case d::syscall_pwrite64:
      return do_write(*state, proc, context, space, /*positional=*/true);
    case d::syscall_close:
      return do_close(*state, proc, context);
    case d::syscall_lseek:
      return do_lseek(proc, context);
    case d::syscall_fstat:
      return do_fstat(*state, proc, context, space);
    case d::syscall_getdents64:
      return do_getdents64(proc, context, space);
    case d::syscall_fchdir:
      return do_fchdir(proc, context);
    case d::syscall_ftruncate:
      return do_ftruncate(proc, context);
    case d::syscall_fchmod:
      return do_fchmod(proc, context);
    case d::syscall_dup:
      return do_dup(*state, proc, context);
    case d::syscall_dup2:
    case d::syscall_dup3:
      return do_dup2(*state, proc, context);
    case d::syscall_fcntl:
      return do_fcntl(*state, proc, context);
    case d::syscall_ioctl:
      return do_ioctl(proc, context);
    case d::syscall_flock:
    case d::syscall_fsync:
    case d::syscall_fdatasync:
      return do_fsync_like(proc, context);
    case d::syscall_fallocate:
      return do_fallocate(proc, context);
    case d::syscall_fstatfs:
      return do_fstatfs(*state, proc, context, space);

    // path-based: authoritative — the memfs is the whole namespace.
    case d::syscall_open:
    case d::syscall_openat:
      return do_open(*state, proc, context, space);
    case d::syscall_stat:
    case d::syscall_lstat:
    case d::syscall_newfstatat:
      return do_stat(*state, proc, context, space);
    case d::syscall_statx:
      return do_statx(*state, proc, context, space);
    case d::syscall_mkdir:
    case d::syscall_mkdirat:
      return do_mkdir(proc, context, space);
    case d::syscall_rmdir:
    case d::syscall_unlink:
    case d::syscall_unlinkat:
      return do_unlink(proc, context, space);
    case d::syscall_rename:
    case d::syscall_renameat:
      return do_rename(proc, context, space);
    case d::syscall_symlink:
    case d::syscall_symlinkat:
      return do_symlink(proc, context, space);
    case d::syscall_link:
    case d::syscall_linkat:
      return do_link(proc, context, space);
    case d::syscall_readlink:
    case d::syscall_readlinkat:
      return do_readlink(proc, context, space);
    case d::syscall_access:
    case d::syscall_faccessat:
    case d::syscall_faccessat2:
      return do_access(proc, context, space);
    case d::syscall_chdir:
      return do_chdir(proc, context, space);
    case d::syscall_getcwd:
      return do_getcwd(proc, context, space);
    case d::syscall_truncate:
      return do_truncate(proc, context, space);
    case d::syscall_chmod:
    case d::syscall_fchmodat:
      return do_chmod_path(proc, context, space);
    case d::syscall_umask:
      return do_umask(proc, context);
    case d::syscall_utime:
      return do_utime(proc, context, space);
    case d::syscall_utimensat:
      return do_utimensat(proc, context, space);
    case d::syscall_statfs:
      return do_statfs(*state, proc, context, space);

    default:
      return std::nullopt;
    }
  } catch (const std::bad_alloc&) {
    return d::return_error(context, d::linux_enomem);
  } catch (...) {
    return d::return_error(context, d::linux_eio);
  }
}

} // namespace x86sim::linux_syscalls
