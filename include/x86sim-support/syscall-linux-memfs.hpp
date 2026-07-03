// Opt-in memfs syscall personality: Linux file syscalls backed by the portable
// in-memory filesystem (x86sim-support/memfs.hpp) instead of the host.
//
// SysMemFs is one aggregate handler satisfying the SyscallHandler concept, so
// it composes with operator| like every other handler. Its dispatch contract:
//   * std::nullopt for syscall numbers it does not cover — normal chain
//     fall-through.
//   * std::nullopt for fd-based syscalls on fds not in its table — no fd
//     number is special; whoever owns that fd (e.g. the Python bindings'
//     stream routing) handles it. Wiring up stdio is the embedder's job.
//   * Path-based syscalls are handled authoritatively once decoded: the
//     sandbox is the whole filesystem view, so a miss is -ENOENT, never a
//     fall-through to the host.
//
// All state is owned by a MemFsState instance shared between the handler and
// the embedder (never TU-global), so successive machines cannot alias each
// other's filesystems.
#ifndef X86SIM_SUPPORT_SYSCALL_LINUX_MEMFS_HPP
#define X86SIM_SUPPORT_SYSCALL_LINUX_MEMFS_HPP

#include "x86sim-support/memfs.hpp"
#include "x86sim-support/namespaces.hpp"
#include "x86sim-support/syscall-linux.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>

namespace x86sim::linux_syscalls {

struct MemFsState : std::enable_shared_from_this<MemFsState> {
  // One open-file description, shared across dup'd fds so the offset is
  // shared per POSIX. `path` is the canonical path at open time, used as the
  // base for dirfd-relative resolution and fchdir.
  struct OpenFile {
    std::shared_ptr<memfs::Node> node;
    std::string path;
    std::uint64_t offset = 0;
    word_t status_flags = 0; // Linux O_* access mode and status flags
    // getdents64 resume key: name of the last returned entry. The two
    // synthetic entries "." and ".." are emitted first, tracked by dir_phase.
    std::string dir_cursor;
    int dir_phase = 0; // 0: next is ".", 1: next is "..", 2: real entries
    std::uint64_t dir_cookie = 0;
  };

  struct FdEntry {
    std::shared_ptr<OpenFile> file;
    bool cloexec = false;
  };

  struct ProcState {
    ProcessNamespaces ns;
    // The fd table starts empty: memfs owns exactly the fds recorded here and
    // nothing else. std::map so allocation scans run in fd order.
    std::map<int, FdEntry> fds;
    // Reserved (embedder-owned) fd numbers the guest has close()d. The number
    // is dead — reads/writes fail with EBADF instead of falling through to the
    // embedder — until a later open/dup reuses it.
    std::set<int> closed_reserved;
    std::string cwd = "/";
    word_t umask = 022;
  };

  explicit MemFsState(std::shared_ptr<memfs::Filesystem> fs);
  ~MemFsState();

  MemFsState(const MemFsState&) = delete;
  MemFsState& operator=(const MemFsState&) = delete;

  // Lazily creates the per-process state (sharing this state's mount
  // namespace) and registers the pid for the file-backed-mmap loader shim.
  [[nodiscard]] ProcState& proc(ProcessId pid);
  // fork hook: child shares namespaces and open-file descriptions, copies the
  // fd table and cwd (matching fork()).
  void clone_proc(ProcessId parent, ProcessId child);
  void erase_proc(ProcessId pid);

  // Lowest fd not used by this table and not reserved by the embedder.
  [[nodiscard]] int allocate_fd(const ProcState& proc, int minimum = 0) const;

  // Host-side API: bind an open memfs file to an arbitrary fd number
  // (including 0/1/2 — this is how an embedder makes stdio a memfs file).
  void install_at(ProcessId pid, int fd, OpenFile file, bool cloexec = false);

  std::shared_ptr<memfs::Filesystem> fs;
  std::shared_ptr<MountNamespace> mount;
  std::unordered_map<ProcessId, ProcState> procs;
  // fd numbers owned by the embedder (e.g. fds routed to Python streams);
  // guest opens never allocate these, but dup2/install_at may shadow them.
  std::set<int> reserved_fds;
};

// Creates the shared state (with a fresh filesystem when none is supplied) and
// installs the file-backed-mmap loader shim on first use.
[[nodiscard]] std::shared_ptr<MemFsState> make_memfs_state(std::shared_ptr<memfs::Filesystem> fs = nullptr);

struct SysMemFs {
  explicit SysMemFs(std::shared_ptr<MemFsState> state) noexcept : state(std::move(state)) {}

  [[nodiscard]] std::optional<SyscallResult> try_syscall(Machine&, ProcessId, CpuState&, AddressSpace&,
                                                         SyscallKind) noexcept;

  std::shared_ptr<MemFsState> state;
};

} // namespace x86sim::linux_syscalls

#endif
