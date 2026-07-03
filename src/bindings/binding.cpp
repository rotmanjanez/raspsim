#include <pybind11/pybind11.h>

#include "x86sim-support/cpuid.hpp"
#include "x86sim-support/memfs.hpp"
#include "x86sim-support/syscall-linux-memfs.hpp"
#include "x86sim-support/syscall-linux.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <format>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace py::literals;

using x86sim::address_t;
using x86sim::Register;
using x86sim::RegisterFile;
using x86sim::word_t;
using x86sim::XmmRegister;
using x86sim::XmmValue;

// Python-facing protection flags. The member names and composites match the old
// binding so the pure-Python layer (elf.py) is untouched; the underlying values
// follow x86sim::Protection semantics.
enum class Prot {
  READ = static_cast<int>(x86sim::Protection::read),
  WRITE = static_cast<int>(x86sim::Protection::write),
  EXEC = static_cast<int>(x86sim::Protection::execute),
  NONE = static_cast<int>(x86sim::Protection::none),
  RW = READ | WRITE,
  RX = READ | EXEC,
  RWX = READ | WRITE | EXEC
};

bool hasProt(Prot p, Prot q) {
  return (static_cast<int>(p) & static_cast<int>(q)) == static_cast<int>(q);
}

Prot addProt(Prot p, Prot q) {
  return static_cast<Prot>(static_cast<int>(p) | static_cast<int>(q));
}

x86sim::Protection toProtection(Prot p) {
  return static_cast<x86sim::Protection>(static_cast<std::uint8_t>(p));
}

// ELF segment flags (PF_R=4, PF_W=2, PF_X=1) -> Prot.
Prot getProtFromELFSegment(int flags) {
  int prot = static_cast<int>(Prot::NONE);
  if (flags & 0x4)
    prot |= static_cast<int>(Prot::READ);
  if (flags & 0x2)
    prot |= static_cast<int>(Prot::WRITE);
  if (flags & 0x1)
    prot |= static_cast<int>(Prot::EXEC);
  return static_cast<Prot>(prot);
}

class PyMachine;

class AddrRef {
public:
  AddrRef(PyMachine* sim, address_t virtaddr) : virtaddr(virtaddr), sim(sim) {}
  AddrRef() : virtaddr(0), sim(nullptr) {}

  operator address_t() const { return virtaddr; }

  AddrRef operator+(address_t offset) { return AddrRef{sim, virtaddr + offset}; }
  AddrRef operator-(address_t offset) { return AddrRef{sim, virtaddr - offset}; }

  bool operator==(const AddrRef& other) const { return virtaddr == other.virtaddr; }
  bool operator!=(const AddrRef& other) const { return virtaddr != other.virtaddr; }
  bool operator<(const AddrRef& other) const { return virtaddr < other.virtaddr; }
  bool operator<=(const AddrRef& other) const { return virtaddr <= other.virtaddr; }
  bool operator>(const AddrRef& other) const { return virtaddr > other.virtaddr; }
  bool operator>=(const AddrRef& other) const { return virtaddr >= other.virtaddr; }

  void write(py::bytes&& bts);
  py::bytes read(word_t size = 1);

  address_t virtaddr;

protected:
  PyMachine* sim;
};

class RaspsimException : public std::exception {
public:
  RaspsimException(const x86sim::X86Exception& exc) : vector(exc.vector), msg(std::format("{}", exc)) {}

  RaspsimException(const char* m) : vector(-1), msg(m) {}
  RaspsimException(std::string m) : vector(-1), msg(std::move(m)) {}

  int getVector() const { return vector; }
  const char* what() const noexcept override { return msg.c_str(); }

private:
  int vector;
  std::string msg;
};

#define RASPSIM_INHERIT_EXCEPTION(name)                                                                                \
  class name : public RaspsimException {                                                                               \
  public:                                                                                                              \
    name(const RaspsimException& e) : RaspsimException(e) {}                                                           \
    name(RaspsimException&& e) : RaspsimException(std::move(e)) {}                                                     \
  }

// X86 exceptions in order of their exception numbers
RASPSIM_INHERIT_EXCEPTION(RaspsimDivideException);
RASPSIM_INHERIT_EXCEPTION(RaspsimDebugException);
RASPSIM_INHERIT_EXCEPTION(RaspsimNMIException);
RASPSIM_INHERIT_EXCEPTION(RaspsimBreakpointException);
RASPSIM_INHERIT_EXCEPTION(RaspsimOverflowException);
RASPSIM_INHERIT_EXCEPTION(RaspsimBoundsException);
RASPSIM_INHERIT_EXCEPTION(RaspsimInvalidOpcodeException);
RASPSIM_INHERIT_EXCEPTION(RaspsimFPUNotAvailException);
RASPSIM_INHERIT_EXCEPTION(RaspsimDoubleFaultException);
RASPSIM_INHERIT_EXCEPTION(RaspsimCoprocOverrunException);
RASPSIM_INHERIT_EXCEPTION(RaspsimInvalidTSSException);
RASPSIM_INHERIT_EXCEPTION(RaspsimSegNotPresentException);
RASPSIM_INHERIT_EXCEPTION(RaspsimStackFaultException);
RASPSIM_INHERIT_EXCEPTION(RaspsimGPFaultException);
RASPSIM_INHERIT_EXCEPTION(RaspsimPageFaultException);
RASPSIM_INHERIT_EXCEPTION(RaspsimSpuriousIntException);
RASPSIM_INHERIT_EXCEPTION(RaspsimFPUException);
RASPSIM_INHERIT_EXCEPTION(RaspsimUnalignedException);
RASPSIM_INHERIT_EXCEPTION(RaspsimMachineCheckException);
RASPSIM_INHERIT_EXCEPTION(RaspsimSSEException);

class RegisterFileRef;

// Configurable dispatcher for the Python binding. By default this is a
// register/memory-level simulator: an int 0x80 signals the guest is done and any
// other syscall is unsupported (surfaced as an exception). Optionally it can:
//   * enable the portable Linux "malloc/free" heap (brk + anonymous
//     mmap/munmap/mremap) when `enable_heap` is set, and
//   * route guest read/write on configured fds to Python file-like objects
//     (e.g. io.BytesIO) instead of any host fd.
// --- memfs host-side view ---------------------------------------------------

// Translate a memfs errno (a Linux value) into a real Python OSError so the
// interpreter picks the right subclass (FileNotFoundError, FileExistsError,
// ...). The common values are identical across POSIX platforms, so subclass
// mapping works for wheels built anywhere.
[[nodiscard]] const char* memfs_error_message(int error) {
  switch (error) {
  case 1:
    return "Operation not permitted";
  case 2:
    return "No such file or directory";
  case 13:
    return "Permission denied";
  case 17:
    return "File exists";
  case 20:
    return "Not a directory";
  case 21:
    return "Is a directory";
  case 22:
    return "Invalid argument";
  case 27:
    return "File too large";
  case 36:
    return "File name too long";
  case 39:
    return "Directory not empty";
  case 40:
    return "Too many levels of symbolic links";
  default:
    return "I/O error";
  }
}

[[noreturn]] void throw_memfs_error(int error, const std::string& path) {
  py::object exc = py::module_::import("builtins").attr("OSError")(error, memfs_error_message(error), path);
  // OSError(errno, ...) normalizes to the matching subclass; raise with the
  // instance's own type so pybind11 sees a consistent (type, value) pair.
  PyErr_SetObject(reinterpret_cast<PyObject*>(Py_TYPE(exc.ptr())), exc.ptr());
  throw py::error_already_set();
}

// Host-side query/prepopulation handle over a memfs. Constructible standalone
// (Fs()) to build a filesystem before any Machine exists, shareable between
// Machines via the Machine(memfs=fs) kwarg, and returned by Machine.fs. Holds
// the shared state, so the tree outlives any Machine using it. Paths are
// resolved from the filesystem root (host-side cwd is always "/").
class PyFs {
public:
  PyFs() : state_(x86sim::linux_syscalls::make_memfs_state()) {}
  explicit PyFs(std::shared_ptr<x86sim::linux_syscalls::MemFsState> state) : state_(std::move(state)) {}

  [[nodiscard]] const std::shared_ptr<x86sim::linux_syscalls::MemFsState>& state() const { return state_; }

  py::list listdir(const std::string& path) {
    auto node = resolve(path);
    if (!node->is_dir())
      throw_memfs_error(20, path); // ENOTDIR
    py::list names;
    for (const auto& [name, child] : node->dir().entries)
      names.append(name);
    return names;
  }

  py::bytes read_bytes(const std::string& path) {
    auto node = resolve(path);
    if (node->is_dir())
      throw_memfs_error(21, path); // EISDIR
    const auto& bytes = node->file().bytes;
    return {reinterpret_cast<const char*>(bytes.data()), bytes.size()};
  }

  void write_bytes(const std::string& path, py::bytes data) {
    auto& fs = *state_->fs;
    std::shared_ptr<x86sim::memfs::Node> node;
    auto resolved = fs.resolve(path, fs.root(), "/");
    if (resolved) {
      node = *resolved;
    } else if (resolved.error() == 2) { // ENOENT: create (parents must exist)
      auto parent = fs.resolve_parent(path, fs.root(), "/");
      if (!parent)
        throw_memfs_error(parent.error(), path);
      auto created = fs.create_file(*parent->dir, parent->leaf, x86sim::memfs::default_file_permissions);
      if (!created)
        throw_memfs_error(created.error(), path);
      node = *created;
    } else {
      throw_memfs_error(resolved.error(), path);
    }
    if (node->is_dir())
      throw_memfs_error(21, path); // EISDIR

    std::string raw = std::move(data);
    if (auto truncated = fs.truncate(*node, 0); !truncated)
      throw_memfs_error(truncated.error(), path);
    if (auto written = fs.write(*node, 0, std::as_bytes(std::span(raw.data(), raw.size()))); !written)
      throw_memfs_error(written.error(), path);
  }

  void mkdir(const std::string& path) {
    auto parent = resolve_parent(path);
    if (auto created = state_->fs->make_directory(*parent.dir, parent.leaf, 0755); !created)
      throw_memfs_error(created.error(), path);
  }

  void unlink(const std::string& path) {
    auto parent = resolve_parent(path);
    if (auto removed = state_->fs->unlink(*parent.dir, parent.leaf); !removed)
      throw_memfs_error(removed.error(), path);
  }

  void rmdir(const std::string& path) {
    auto parent = resolve_parent(path);
    if (auto removed = state_->fs->remove_directory(*parent.dir, parent.leaf); !removed)
      throw_memfs_error(removed.error(), path);
  }

  void rename(const std::string& source, const std::string& target) {
    auto source_parent = resolve_parent(source);
    auto target_parent = resolve_parent(target);
    if (auto renamed = state_->fs->rename(*source_parent.dir, source_parent.leaf, *target_parent.dir,
                                          target_parent.leaf);
        !renamed)
      throw_memfs_error(renamed.error(), source);
  }

  void symlink(const std::string& target, const std::string& link_path) {
    auto parent = resolve_parent(link_path);
    if (auto created = state_->fs->make_symlink(*parent.dir, parent.leaf, target); !created)
      throw_memfs_error(created.error(), link_path);
  }

  std::string readlink(const std::string& path) {
    auto node = resolve(path, /*follow=*/false);
    if (!node->is_symlink())
      throw_memfs_error(22, path); // EINVAL
    return node->symlink().target;
  }

  py::dict stat(const std::string& path, bool follow_symlinks) {
    auto node = resolve(path, follow_symlinks);
    const auto info = state_->fs->stat(*node);
    py::dict result;
    result["st_mode"] = info.mode;
    result["st_ino"] = info.ino;
    result["st_nlink"] = info.nlink;
    result["st_size"] = info.size;
    result["st_atime"] = info.atime;
    result["st_mtime"] = info.mtime;
    result["st_ctime"] = info.ctime;
    return result;
  }

  void truncate(const std::string& path, std::uint64_t size) {
    auto node = resolve(path);
    if (auto truncated = state_->fs->truncate(*node, size); !truncated)
      throw_memfs_error(truncated.error(), path);
  }

  bool exists(const std::string& path) { return try_resolve(path, true) != nullptr; }
  bool is_file(const std::string& path) {
    auto node = try_resolve(path, true);
    return node && node->is_file();
  }
  bool is_dir(const std::string& path) {
    auto node = try_resolve(path, true);
    return node && node->is_dir();
  }
  bool is_symlink(const std::string& path) {
    auto node = try_resolve(path, false);
    return node && node->is_symlink();
  }

private:
  [[nodiscard]] std::shared_ptr<x86sim::memfs::Node> resolve(const std::string& path, bool follow = true) {
    auto& fs = *state_->fs;
    auto resolved = fs.resolve(path, fs.root(), "/", follow);
    if (!resolved)
      throw_memfs_error(resolved.error(), path);
    return *resolved;
  }

  [[nodiscard]] std::shared_ptr<x86sim::memfs::Node> try_resolve(const std::string& path, bool follow) {
    auto& fs = *state_->fs;
    auto resolved = fs.resolve(path, fs.root(), "/", follow);
    return resolved ? *resolved : nullptr;
  }

  [[nodiscard]] x86sim::memfs::Filesystem::ParentRef resolve_parent(const std::string& path) {
    auto& fs = *state_->fs;
    auto parent = fs.resolve_parent(path, fs.root(), "/");
    if (!parent)
      throw_memfs_error(parent.error(), path);
    return *parent;
  }

  std::shared_ptr<x86sim::linux_syscalls::MemFsState> state_;
};

// Host I/O is never used: stdin/stdout/stderr map to caller-supplied Python
// objects only.
class PyHost : public x86sim::HostCallbacks {
public:
  // The portable heap handlers are constructed fresh per PyHost (hence per
  // PyMachine) so SysMmap::next_mapping_address resets for each Machine. The pid
  // is unique per PyHost so the static brk_states map in the support library
  // never reuses heap state across successive Machine() instances.
  PyHost() : pid(next_pid()) {}

  x86sim::SyscallResult syscall(x86sim::Machine& machine, x86sim::CpuState& ctx, x86sim::AddressSpace& space,
                                x86sim::SyscallKind kind) override {
    using x86sim::StopReason;
    namespace abi = x86sim::linux_syscalls::abi;

    // int 0x80 stays the guest-exit sentinel: existing README examples rely on
    // `mov rax, -1; int 0x80` to stop the simulator.
    if (kind == x86sim::SyscallKind::int80)
      return guest_exit();

    // The portable handlers (and the ABI helpers) only decode the 64-bit
    // `syscall` instruction; anything else is unsupported.
    if (kind != x86sim::SyscallKind::syscall64)
      return unsupported();

    const word_t n = abi::syscall_number(ctx);

    if (n == abi::exit || n == abi::exit_group)
      return guest_exit();

    // Opt-in memfs personality, consulted before the stream shortcuts below so
    // a memfs entry installed at any fd number (guest dup2 or host map_fd)
    // shadows a Python stream at the same number. Fds that are not in the
    // memfs table fall through here (no fd number is special to the memfs).
    if (memfs_handler) {
      if (auto r = memfs_handler->try_syscall(machine, pid, ctx, space, kind))
        return *r;
      // Metadata syscalls on the embedder's own stream fds are answered here:
      // glibc stdio startup fstat()s and ioctl()s fd 1, and nothing else in
      // the wheel implements those.
      if (n == kSyscallFstat) {
        if (auto r = stream_fstat(ctx, space))
          return *r;
      }
      if (n == kSyscallIoctl) {
        if (auto r = stream_ioctl(ctx))
          return *r;
      }
    }

    if (n == abi::read)
      return do_read(ctx, space);
    if (n == abi::write)
      return do_write(ctx, space);
    if (n == kSyscallReadlink || n == kSyscallReadlinkat)
      return do_readlink(ctx, space);

    // Opt-in portable Linux glibc-startup chain. Every handler self-decodes the
    // syscall number and returns std::nullopt for anything it does not handle,
    // so consulting the chain is harmless for unrelated syscalls.
    if (enable_glibc) {
      if (auto r = glibc_syscalls.try_syscall(machine, pid, ctx, space, kind))
        return *r;
    }

    return unsupported();
  }

  x86sim::CpuidResult cpuid(x86sim::Machine&, x86sim::CpuState&, x86sim::AddressSpace&,
                            x86sim::CpuidRequest request) noexcept override {
    return x86sim::defaults::default_cpuid(request);
  }

  // If an I/O callback raised during the last run, re-raise the original Python
  // exception (clearing the stash). Called by PyMachine::run after run() returns.
  void rethrow_pending_error() {
    if (pending_error)
      std::rethrow_exception(std::exchange(pending_error, nullptr));
  }

  // Enables the portable Linux glibc-startup chain: the malloc/free heap
  // (brk + anonymous mmap/munmap/mremap) plus arch_prctl, set_tid_address,
  // set_robust_list, rseq, prlimit64, uname, futex and the synthetic signal
  // syscalls. All of these are register/memory-only handlers from the support
  // library (no host syscalls).
  bool enable_glibc = false;
  // Guest fd -> Python file-like object. Seeded with 0/1/2 in the ctor; an
  // absent (or None) entry means that fd is unconfigured.
  std::unordered_map<int, py::object> fds;
  // Optional Python callable resolving a symlink path to its target,
  // readlink(path: str) -> str | None. glibc reads /proc/self/exe at startup.
  // None ⇒ readlink/readlinkat return -ENOSYS; a callback returning None ⇒
  // -ENOENT. Routed through Python (like read/write) rather than the host so the
  // binding stays free of real filesystem access.
  py::object readlink_cb = py::none();
  // Opt-in memfs personality: null/empty when disabled. The state is shared
  // with the PyFs handle(s) the user sees, so the tree outlives the Machine.
  std::shared_ptr<x86sim::linux_syscalls::MemFsState> memfs_state;
  std::optional<x86sim::linux_syscalls::SysMemFs> memfs_handler;

  [[nodiscard]] x86sim::linux_syscalls::ProcessId process_id() const { return pid; }

private:
  // Cap a single read/write transfer so a bogus guest count cannot trigger a
  // huge allocation; oversized requests become a (legal) short transfer. This
  // mirrors the spirit of the bounds checks in the syscall support library.
  static constexpr word_t kMaxTransfer = word_t{1} << 31; // 2 GiB

  // x86-64 syscall numbers handled directly here (routed to Python), and the
  // negated Linux errno values returned to the guest on the failure paths.
  static constexpr word_t kSyscallReadlink = 89;
  static constexpr word_t kSyscallReadlinkat = 267;
  static constexpr word_t kSyscallFstat = 5;
  static constexpr word_t kSyscallIoctl = 16;
  static constexpr std::int64_t kLinuxEnoent = 2;
  static constexpr std::int64_t kLinuxEnosys = 38;
  static constexpr std::int64_t kLinuxEnotty = 25;

  static x86sim::SyscallResult guest_exit() {
    return {.reason = x86sim::StopReason::guest_exit, .continue_execution = false, .message = {}};
  }

  static x86sim::SyscallResult unsupported() {
    return {.reason = x86sim::StopReason::unsupported_syscall,
            .continue_execution = false,
            .message = "Syscall not supported"};
  }

  static x86sim::linux_syscalls::ProcessId next_pid() {
    static std::atomic<x86sim::linux_syscalls::ProcessId> counter{x86sim::linux_syscalls::initial_process_id};
    return counter.fetch_add(1, std::memory_order_relaxed);
  }

  // Stop the run and stash a Python exception raised by an I/O callback so
  // PyMachine::run can re-raise the original Python exception to the caller
  // instead of letting it unwind through the simulator core mid-syscall.
  x86sim::SyscallResult py_error_stop() {
    pending_error = std::current_exception();
    return {x86sim::StopReason::host_request, false, "Python I/O callback raised an exception"};
  }

  // The GIL is held throughout Machine::run (pybind11 holds it by default and we
  // never release it), so invoking the Python .read/.write callbacks below is
  // safe without any additional acquire. A callback that raises is captured by
  // py_error_stop() and re-raised after run() returns.
  x86sim::SyscallResult do_read(x86sim::CpuState& ctx, x86sim::AddressSpace& space) {
    namespace abi = x86sim::linux_syscalls::abi;
    const int fd = static_cast<int>(abi::syscall_arg(ctx, 0));
    auto it = fds.find(fd);
    if (it == fds.end() || it->second.is_none() || !py::hasattr(it->second, "read"))
      return unsupported();

    const address_t buffer = abi::syscall_arg(ctx, 1);
    const word_t count = std::min<word_t>(abi::syscall_arg(ctx, 2), kMaxTransfer);

    std::string data;
    try {
      py::object result = it->second.attr("read")(static_cast<std::size_t>(count));
      // Accept any bytes-like return (bytes/bytearray/memoryview, ...).
      PyObject* raw_bytes = PyBytes_FromObject(result.ptr());
      if (!raw_bytes)
        throw py::error_already_set();
      data = static_cast<std::string>(py::reinterpret_steal<py::bytes>(raw_bytes));
    } catch (py::error_already_set&) {
      return py_error_stop();
    }

    // Honour the configured count (handles short reads and over-long returns).
    const std::size_t n = std::min<std::size_t>(data.size(), static_cast<std::size_t>(count));
    if (n > 0) {
      auto bytes = std::as_bytes(std::span(data.data(), n));
      if (auto r = space.write(buffer, bytes); !r)
        return unsupported();
    }
    return abi::return_value(ctx, static_cast<std::int64_t>(n));
  }

  x86sim::SyscallResult do_write(x86sim::CpuState& ctx, x86sim::AddressSpace& space) {
    namespace abi = x86sim::linux_syscalls::abi;
    const int fd = static_cast<int>(abi::syscall_arg(ctx, 0));
    auto it = fds.find(fd);
    if (it == fds.end() || it->second.is_none() || !py::hasattr(it->second, "write"))
      return unsupported();

    const address_t buffer = abi::syscall_arg(ctx, 1);
    const word_t count = std::min<word_t>(abi::syscall_arg(ctx, 2), kMaxTransfer);

    // Stream the guest buffer to the Python object in bounded chunks: this caps
    // our own allocation regardless of `count`, and honours short writes (a
    // file-like .write() may consume fewer bytes than offered) so the value
    // returned to the guest matches the write(2) ABI.
    py::object writer = it->second.attr("write");
    std::array<std::byte, kIoChunk> chunk{};
    word_t total = 0;
    try {
      while (total < count) {
        const std::size_t want = static_cast<std::size_t>(std::min<word_t>(count - total, chunk.size()));
        if (auto r = space.read(buffer + total, std::span(chunk.data(), want)); !r)
          return total == 0 ? unsupported() : abi::return_value(ctx, static_cast<std::int64_t>(total));

        py::object ret = writer(py::bytes(reinterpret_cast<const char*>(chunk.data()), want));
        // A binary stream returns the number of bytes consumed; None means "all"
        // (matching e.g. text/raw stream conventions). Clamp defensively.
        std::size_t written = want;
        if (py::isinstance<py::int_>(ret)) {
          const auto reported = ret.cast<long long>();
          written = reported <= 0 ? 0 : std::min<std::size_t>(static_cast<std::size_t>(reported), want);
        }
        total += static_cast<word_t>(written);
        if (written < want) // short write: stop and report the partial count
          break;
      }
    } catch (py::error_already_set&) {
      return total == 0 ? py_error_stop() : abi::return_value(ctx, static_cast<std::int64_t>(total));
    }
    return abi::return_value(ctx, static_cast<std::int64_t>(total));
  }

  // readlink(path, buf, size) / readlinkat(dirfd, path, buf, size). Resolution
  // is delegated to the optional Python `readlink_cb` so the binding performs no
  // real filesystem access; the dirfd of readlinkat is ignored (paths are
  // resolved by the callback). Unlike read(2), readlink(2) does not append a NUL
  // and the result is truncated (not an error) when it exceeds the buffer.
  x86sim::SyscallResult do_readlink(x86sim::CpuState& ctx, x86sim::AddressSpace& space) {
    namespace abi = x86sim::linux_syscalls::abi;
    if (readlink_cb.is_none())
      return abi::return_value(ctx, -kLinuxEnosys);

    const bool is_at = abi::syscall_number(ctx) == kSyscallReadlinkat;
    const address_t path_addr = abi::syscall_arg(ctx, is_at ? 1 : 0);
    const address_t buffer = abi::syscall_arg(ctx, is_at ? 2 : 1);
    const word_t size = std::min<word_t>(abi::syscall_arg(ctx, is_at ? 3 : 2), kMaxTransfer);

    std::string path;
    if (!read_guest_cstring(space, path_addr, path))
      return unsupported();

    std::string target;
    try {
      py::object result = readlink_cb(path);
      if (result.is_none())
        return abi::return_value(ctx, -kLinuxEnoent);
      target = result.cast<std::string>();
    } catch (py::error_already_set&) {
      return py_error_stop();
    }

    const std::size_t n = std::min<std::size_t>(target.size(), static_cast<std::size_t>(size));
    if (n > 0) {
      auto bytes = std::as_bytes(std::span(target.data(), n));
      if (auto r = space.write(buffer, bytes); !r)
        return unsupported();
    }
    return abi::return_value(ctx, static_cast<std::int64_t>(n));
  }

  // fstat on one of the embedder's stream fds: answer with a synthetic
  // character-device stat (mode 020620, like a tty) so glibc stdio startup
  // succeeds. Only used when memfs is enabled; returns std::nullopt for fds
  // that are not configured streams.
  std::optional<x86sim::SyscallResult> stream_fstat(x86sim::CpuState& ctx, x86sim::AddressSpace& space) {
    namespace abi = x86sim::linux_syscalls::abi;
    const int fd = static_cast<int>(abi::syscall_arg(ctx, 0));
    auto it = fds.find(fd);
    if (it == fds.end() || it->second.is_none())
      return std::nullopt;

    std::array<std::byte, 144> bytes{};
    auto write_le = [&](std::size_t offset, std::uint64_t value, std::size_t width) noexcept {
      for (std::size_t i = 0; i < width; ++i)
        bytes[offset + i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);
    };
    write_le(0, 1, 8);       // st_dev
    write_le(8, 1, 8);       // st_ino
    write_le(16, 1, 8);      // st_nlink
    write_le(24, 020620, 4); // st_mode: S_IFCHR | 0620
    write_le(56, 1024, 8);   // st_blksize

    if (auto r = space.write(abi::syscall_arg(ctx, 1), bytes); !r)
      return unsupported();
    return abi::return_value(ctx, 0);
  }

  std::optional<x86sim::SyscallResult> stream_ioctl(x86sim::CpuState& ctx) {
    namespace abi = x86sim::linux_syscalls::abi;
    const int fd = static_cast<int>(abi::syscall_arg(ctx, 0));
    auto it = fds.find(fd);
    if (it == fds.end() || it->second.is_none())
      return std::nullopt;
    return abi::return_value(ctx, -kLinuxEnotty);
  }

  // Read a NUL-terminated guest string into `out` (without the terminator).
  // Returns false on unreadable memory or a missing terminator within the cap.
  static bool read_guest_cstring(x86sim::AddressSpace& space, address_t addr, std::string& out) {
    constexpr std::size_t kMaxPath = 4096;
    out.clear();
    for (std::size_t i = 0; i < kMaxPath; ++i) {
      std::byte byte{};
      if (auto r = space.read(addr + i, std::span(&byte, 1)); !r)
        return false;
      char c = static_cast<char>(byte);
      if (c == '\0')
        return true;
      out.push_back(c);
    }
    return false;
  }

  // Per-chunk transfer size for streamed writes; bounds our own buffering.
  static constexpr std::size_t kIoChunk = 64 * 1024;

  std::exception_ptr pending_error;
  x86sim::linux_syscalls::ProcessId pid;
  // The portable glibc-startup syscalls. Built fresh per PyHost so the heap
  // handlers' next_mapping_address / break state reset for each Machine. read,
  // write, exit/exit_group and readlink are handled above (routed to Python), so
  // they are intentionally absent from this chain. SysGetIdentity is included
  // (getpid/getuid/...; glibc probes uid/gid at startup when AT_SECURE is
  // absent); it returns synthetic constants and the pid passed below, and does
  // not touch the (unused) ProcessTable. No other ProcessTable-backed handlers
  // (fork/wait/...) are wired in.
  decltype(x86sim::linux_syscalls::SysBrk{} | x86sim::linux_syscalls::SysMmap{} | x86sim::linux_syscalls::SysMunmap{} |
           x86sim::linux_syscalls::SysMprotect{} | x86sim::linux_syscalls::SysMremap{} |
           x86sim::linux_syscalls::SysArchPrctl{} | x86sim::linux_syscalls::SysSetTidAddress{} |
           x86sim::linux_syscalls::SysSetRobustList{} | x86sim::linux_syscalls::SysRseq{} |
           x86sim::linux_syscalls::SysPrlimit64{} | x86sim::linux_syscalls::SysUname{} |
           x86sim::linux_syscalls::SysGetIdentity{} | x86sim::linux_syscalls::SysFutex{} |
           x86sim::linux_syscalls::SysSignals{} | x86sim::linux_syscalls::SysGetrandom{}) glibc_syscalls =
      x86sim::linux_syscalls::SysBrk{} | x86sim::linux_syscalls::SysMmap{} | x86sim::linux_syscalls::SysMunmap{} |
      x86sim::linux_syscalls::SysMprotect{} | x86sim::linux_syscalls::SysMremap{} |
      x86sim::linux_syscalls::SysArchPrctl{} | x86sim::linux_syscalls::SysSetTidAddress{} |
      x86sim::linux_syscalls::SysSetRobustList{} | x86sim::linux_syscalls::SysRseq{} |
      x86sim::linux_syscalls::SysPrlimit64{} | x86sim::linux_syscalls::SysUname{} |
      x86sim::linux_syscalls::SysGetIdentity{} | x86sim::linux_syscalls::SysFutex{} |
      x86sim::linux_syscalls::SysSignals{} | x86sim::linux_syscalls::SysGetrandom{};
};

class PyMachine {
public:
  // Build the machine options up front so the Machine (and the AddressSpace
  // bound to it) can be constructed in the member initializer list.
  static x86sim::Options build_options(const char* logfile, bool sse, bool x87, bool perfect_cache,
                                       bool static_branchpred, const char* core) {
    x86sim::Options options;
    options.sse = sse;
    options.x87 = x87;
    options.debug.perfect_cache = perfect_cache;
    options.debug.static_branchpred = static_branchpred;
    options.log.log_filename = logfile;

    // Core model selection. The out-of-order core is the default; the sequential
    // core is slower but executes unaligned loads/stores correctly, which the
    // out-of-order core's alignment-fixup path currently cannot (it stalls until
    // the deadlock detector aborts the run). Code that relies on unaligned SSE
    // accesses -- notably glibc's string routines -- needs core="seq".
    const std::string_view core_name(core);
    if (core_name == "seq" || core_name == "sequential")
      options.core = x86sim::CoreModel::sequential;
    else if (core_name == "ooo" || core_name == "out_of_order")
      options.core = x86sim::CoreModel::out_of_order;
    else
      throw py::value_error(R"(core must be one of "ooo"/"out_of_order" or "seq"/"sequential")");

    return options;
  }

  PyMachine(const char* logfile, bool sse, bool x87, bool perfect_cache, bool static_branchpred, bool glibc,
            py::object stdin_obj, py::object stdout_obj, py::object stdout_err, py::object readlink_cb,
            const char* core, py::object memfs)
      : machine(std::make_unique<x86sim::Machine>(
            host, build_options(logfile, sse, x87, perfect_cache, static_branchpred, core))),
        address_space(*machine) {
    host.enable_glibc = glibc;

    // Opt-in memfs: True creates a fresh filesystem, an Fs object shares an
    // existing one (prepopulation / cross-Machine sharing). memfs owns the
    // whole path namespace, so the readlink callback is mutually exclusive.
    const bool memfs_off = memfs.is_none() || (py::isinstance<py::bool_>(memfs) && !memfs.cast<bool>());
    if (!memfs_off) {
      if (!readlink_cb.is_none())
        throw py::value_error("memfs and readlink are mutually exclusive: with memfs enabled the in-memory "
                              "filesystem resolves readlink()");
      if (py::isinstance<py::bool_>(memfs))
        host.memfs_state = x86sim::linux_syscalls::make_memfs_state();
      else if (py::isinstance<PyFs>(memfs))
        host.memfs_state = memfs.cast<PyFs&>().state();
      else
        throw py::value_error("memfs must be a bool or an x86sim Fs object");
      host.memfs_handler.emplace(host.memfs_state);
    }

    if (!readlink_cb.is_none()) {
      if (!py::hasattr(readlink_cb, "__call__"))
        throw py::value_error("readlink must be a callable: readlink(path: str) -> str | None");
      host.readlink_cb = std::move(readlink_cb);
    }
    // Map guest fds 0/1/2 to the supplied Python file-like objects (None leaves
    // the fd unconfigured, so the guest sees an unsupported syscall on it).
    map_stream(0, std::move(stdin_obj), "read", "stdin");
    map_stream(1, std::move(stdout_obj), "write", "stdout");
    map_stream(2, std::move(stdout_err), "write", "stderr");
    // Every fd routed to a Python stream is reserved: guest opens never
    // allocate those numbers (though dup2/map_fd may deliberately shadow them).
    if (host.memfs_state)
      for (const auto& [fd, stream] : host.fds)
        host.memfs_state->reserved_fds.insert(fd);
  }

  ~PyMachine() = default;

  static constexpr word_t getPageSize() { return x86sim::AddressSpace::kPageSize; }

  x86sim::Machine& m() { return *machine; }
  // The architectural state and address space are caller-owned (one each here,
  // modelling a single guest process) and handed to Machine::run.
  x86sim::CpuState& state() { return cpu_state; }
  x86sim::AddressSpace& space() { return address_space; }

  RegisterFileRef getRegisters();

  AddrRef memmap(address_t start, Prot prot, word_t length, py::bytes data);

  std::size_t cycles() { return machine->stats().cycles; }
  std::size_t instructions() { return machine->stats().instructions; }

  std::string str() { return std::format("{}", static_cast<x86sim::RegisterFile&>(cpu_state)); }

  void run(unsigned long long ninstr);

  // Host-side handle to the in-memory filesystem (only with memfs enabled).
  PyFs fs() {
    if (!host.memfs_state)
      throw py::value_error("memfs is not enabled on this Machine (construct it with memfs=True)");
    return PyFs(host.memfs_state);
  }

  // Bind a memfs file to an arbitrary guest fd number — including 0/1/2, which
  // is how stdio becomes a memfs file. mode follows open(): "r", "w", "a",
  // optionally with "+" (and an ignored "b").
  int map_fd(int fd, const std::string& path, const std::string& mode) {
    if (!host.memfs_state)
      throw py::value_error("memfs is not enabled on this Machine (construct it with memfs=True)");
    if (fd < 0)
      throw py::value_error("fd must be non-negative");

    bool readable = false, writable = false, append = false, create = false, truncate = false;
    int primaries = 0;
    for (const char c : mode) {
      switch (c) {
      case 'r':
        readable = true;
        ++primaries;
        break;
      case 'w':
        writable = create = truncate = true;
        ++primaries;
        break;
      case 'a':
        writable = create = append = true;
        ++primaries;
        break;
      case '+':
        readable = writable = true;
        break;
      case 'b':
      case 't':
        break;
      default:
        throw py::value_error(std::format("invalid mode: '{}'", mode));
      }
    }
    if (primaries != 1)
      throw py::value_error(std::format("invalid mode: '{}'", mode));

    auto& state = *host.memfs_state;
    auto& fs = *state.fs;
    std::shared_ptr<x86sim::memfs::Node> node;
    auto resolved = fs.resolve(path, fs.root(), "/");
    if (resolved) {
      node = *resolved;
    } else if (resolved.error() == 2 && create) { // ENOENT
      auto parent = fs.resolve_parent(path, fs.root(), "/");
      if (!parent)
        throw_memfs_error(parent.error(), path);
      auto created = fs.create_file(*parent->dir, parent->leaf, x86sim::memfs::default_file_permissions);
      if (!created)
        throw_memfs_error(created.error(), path);
      node = *created;
    } else {
      throw_memfs_error(resolved.error(), path);
    }
    if (node->is_dir())
      throw_memfs_error(21, path); // EISDIR: map_fd is for regular files
    if (truncate)
      if (auto truncated = fs.truncate(*node, 0); !truncated)
        throw_memfs_error(truncated.error(), path);

    // Linux O_* encoding: O_WRONLY=1, O_RDWR=2, O_APPEND=02000.
    word_t flags = readable && writable ? 2 : writable ? 1 : 0;
    if (append)
      flags |= 02000;
    state.install_at(host.process_id(), fd,
                     x86sim::linux_syscalls::MemFsState::OpenFile{.node = std::move(node), .status_flags = flags});
    return fd;
  }

  PyHost host;
  x86sim::CpuState cpu_state;
  // Declared (and therefore constructed) before address_space, which binds to it.
  std::unique_ptr<x86sim::Machine> machine;
  x86sim::AddressSpace address_space;

private:
  // Validate that a supplied stream exposes the attribute the guest will need,
  // then register it for the given fd. A None object is left unconfigured.
  void map_stream(int fd, py::object obj, const char* attr, const char* name) {
    if (obj.is_none())
      return;
    if (!py::hasattr(obj, attr))
      throw py::value_error(std::format("{} must be a file-like object exposing a .{}() method", name, attr));
    host.fds[fd] = std::move(obj);
  }
};

void AddrRef::write(py::bytes&& bts) {
  std::string mem{std::move(bts)};
  auto bytes = std::as_bytes(std::span(mem.data(), mem.size()));
  if (auto r = sim->space().write(virtaddr, bytes); !r)
    throw py::value_error(std::format("Trying to write to unmapped memory at {:x}: {}", virtaddr, r.error()));
}

py::bytes AddrRef::read(word_t size) {
  std::string out(size, '\0');
  auto bytes = std::as_writable_bytes(std::span(out.data(), out.size()));
  if (auto r = sim->space().read(virtaddr, bytes); !r)
    throw py::value_error(std::format("Trying to read from unmapped memory at {:x}: {}", virtaddr, r.error()));
  return {std::move(out)};
}

AddrRef PyMachine::memmap(address_t start, Prot prot, word_t length, py::bytes data) {
  std::string bytes = std::move(data);
  if (bytes.size() > 0) {
    if (length > 0 && bytes.size() > length)
      throw py::value_error("Data size must be less than or equal to length");
    length = std::max(static_cast<word_t>(length), static_cast<word_t>(bytes.size()));
  }

  if (length == 0)
    throw py::value_error("Cannot map zero length data");

  word_t offset = start % getPageSize();
  if (auto r = address_space.map(start - offset, length + offset, toProtection(prot)); !r)
    throw py::value_error(std::format("Cannot map memory at {:x}: {}", start, r.error()));

  if (!bytes.empty()) {
    auto raw = std::as_bytes(std::span(bytes.data(), bytes.size()));
    if (auto r = address_space.write(start, raw); !r)
      throw py::value_error(std::format("Cannot write mapped memory at {:x}: {}", start, r.error()));
  }
  return AddrRef(this, start);
}

#define THROW(i, cls)                                                                                                  \
  case i:                                                                                                              \
    throw cls(std::move(exc));

void PyMachine::run(unsigned long long ninstr) {
  x86sim::RunOptions run_options;
  if (ninstr != static_cast<unsigned long long>(-1))
    run_options.instruction_limit = ninstr;

  x86sim::RunResult result = machine->run(cpu_state, address_space, run_options);

  // If a Python I/O callback raised, surface the original Python exception
  // rather than the generic host_request stop below.
  host.rethrow_pending_error();

  using x86sim::StopReason;
  switch (result.reason) {
  case StopReason::guest_exit:
    return;
  case StopReason::instruction_limit:
    throw py::stop_iteration("Reached instruction limit");
  case StopReason::x86_exception: {
    RaspsimException exc =
        result.x86_exception ? RaspsimException(*result.x86_exception) : RaspsimException("Unknown x86 exception");
    switch (result.x86_exception ? static_cast<int>(result.x86_exception->vector) : -1) {
      THROW(0, RaspsimDivideException)
      THROW(1, RaspsimDebugException)
      THROW(2, RaspsimNMIException)
      THROW(3, RaspsimBreakpointException)
      THROW(4, RaspsimOverflowException)
      THROW(5, RaspsimBoundsException)
      THROW(6, RaspsimInvalidOpcodeException)
      THROW(7, RaspsimFPUNotAvailException)
      THROW(8, RaspsimDoubleFaultException)
      THROW(9, RaspsimCoprocOverrunException)
      THROW(10, RaspsimInvalidTSSException)
      THROW(11, RaspsimSegNotPresentException)
      THROW(12, RaspsimStackFaultException)
      THROW(13, RaspsimGPFaultException)
      THROW(14, RaspsimPageFaultException)
      THROW(15, RaspsimSpuriousIntException)
      THROW(16, RaspsimFPUException)
      THROW(17, RaspsimUnalignedException)
      THROW(18, RaspsimMachineCheckException)
      THROW(19, RaspsimSSEException)
    default:
      throw exc;
    }
  }
  case StopReason::unsupported_syscall:
  case StopReason::host_request:
    throw RaspsimException(result.message.empty() ? "Simulation stopped" : result.message);
  }
}

// Map a register name to its enum, mirroring src/raspsim/raspsim.cpp.
std::optional<Register> registerFromName(std::string_view name) {
  using enum Register;
  static constexpr std::pair<std::string_view, Register> names[] = {
      {"rax", rax}, {"rcx", rcx}, {"rdx", rdx}, {"rbx", rbx}, {"rsp", rsp}, {"rbp", rbp},
      {"rsi", rsi}, {"rdi", rdi}, {"r8", r8},   {"r9", r9},   {"r10", r10}, {"r11", r11},
      {"r12", r12}, {"r13", r13}, {"r14", r14}, {"r15", r15}, {"rip", rip}, {"flags", flags},
  };
  for (auto [reg_name, reg] : names) {
    if (name == reg_name)
      return reg;
  }
  return std::nullopt;
}

class XMMRegister {
public:
  XMMRegister(PyMachine* sim, XmmRegister reg) : sim(sim), reg(reg) {}

  template<typename T, typename... Ts>
  std::tuple<T, Ts...> getPacked() {
    constexpr int NItems = sizeof...(Ts) + 1;
    static_assert(sizeof(T) <= 8, "XMM Register can only be cast to vector types with an "
                                  "element type size less than or equal to 64 bits");
    static_assert((NItems * sizeof(T) == 16) || (NItems == 1 && sizeof(T) <= 8),
                  "XMM Register can only be cast to vector types with size equal to "
                  "128 bits or scalar types with a size less than or equal to 64 "
                  "bits");

    std::tuple<T, Ts...> result;
    XmmValue value = sim->state()[reg];
    std::array<word_t, 2> raw = {value.lo, value.hi};
    std::memcpy(static_cast<void*>(&std::get<0>(result)), raw.data(), NItems * sizeof(T));
    return result;
  }

  template<typename T>
  T getSingle() {
    return std::get<0>(getPacked<T>());
  }

  template<typename T, typename... Ts>
  void setPacked(std::tuple<T, Ts...> value) {
    constexpr int NItems = sizeof...(Ts) + 1;
    static_assert(sizeof(T) <= 8, "XMM Register can only be cast to vector types with an "
                                  "element type size less than or equal to 64 bits");
    static_assert((NItems * sizeof(T) == 16) || (NItems == 1 && sizeof(T) <= 8),
                  "XMM Register can only be cast to vector types with size equal to "
                  "128 bits or scalar types with a size less than or equal to 64 "
                  "bits");

    std::array<word_t, 2> raw = {0, 0};
    std::memcpy(raw.data(), static_cast<void*>(&std::get<0>(value)), NItems * sizeof(T));
    sim->state()[reg] = XmmValue{raw[0], raw[1]};
  }

  template<typename T>
  void setSingle(T value) {
    setPacked<T>(std::tuple<T>{value});
  }

private:
  PyMachine* sim;
  XmmRegister reg;
};

class MemImg {
public:
  MemImg(PyMachine& sim) : sim(&sim) {}

  py::bytes getitem(const py::slice& s) const {
    auto startstop = validateSlice(s);
    return AddrRef(sim, startstop.first).read(startstop.second - startstop.first);
  }

  void setitem(const py::slice& s, py::bytes data) {
    auto startstop = validateSlice(s);
    std::string raw{std::move(data)};
    if (raw.size() != startstop.second - startstop.first)
      throw py::value_error("Data size must match slice size");
    AddrRef(sim, startstop.first).write(py::bytes(raw));
  }

private:
  std::pair<address_t, word_t> validateSlice(const py::slice& s) const {
    if (!py::isinstance<py::none>(s.attr("step")))
      throw py::value_error("Step is not supported");
    return {s.attr("start").cast<address_t>(), s.attr("stop").cast<word_t>()};
  }

  PyMachine* sim;
};

class RegisterFileRef {
  friend class XMMRegister;

public:
  RegisterFileRef() = delete;
  RegisterFileRef(PyMachine* sim) : sim(sim) {}

  word_t getRegister(std::string&& regname) {
    auto reg = registerFromName(regname);
    if (!reg)
      throw py::value_error(std::string("Invalid register name '") + regname + "'");
    return sim->state()[*reg];
  }

  void setRegister(std::string&& regname, word_t value) {
    auto reg = registerFromName(regname);
    if (!reg)
      throw py::value_error(std::string("Invalid register name '") + regname + "'");
    sim->state()[*reg] = value;
  }

  template<typename T, int Bits, int Offset = 0>
  inline T getGPRegImpl(Register reg) {
    static_assert(Bits <= 64, "Register size must be less than or equal to 64 bits");
    static_assert(Offset % 8 == 0, "Offset must be a multiple of 8");
    static_assert(Bits % 8 == 0, "Bits must be a multiple of 8");

    word_t mask = Bits == 64 ? ~word_t{0} : ((word_t{1} << Bits) - 1);
    return (sim->state()[reg] >> Offset) & mask;
  }

  template<typename T, unsigned Bits, unsigned Offset = 0, bool ZeroHigh = false>
  inline void setGPRegImpl(Register reg, T value) {
    static_assert(Bits <= 64, "Register size must be less than or equal to 64 bits");
    static_assert(Offset + Bits <= 64, "Offset + Bits must be less than or equal to 64");
    static_assert(Offset % 8 == 0, "Offset must be a multiple of 8");
    static_assert(Bits % 8 == 0, "Bits must be a multiple of 8");
    static_assert(!(ZeroHigh && Bits != 32 && Offset != 0),
                  "ZeroHigh can only be used with 32-bit registers at offset 0");
    static_assert(!(Offset > 0 && Bits != 8), "Offset can only be used with 8-bit registers");

    word_t current = sim->state()[reg];
    int rshift = ZeroHigh ? 0 : 64 - Bits;

    // mask contains 1s in the bits that are not part of the register that
    // should be modified
    word_t mask = ((~word_t{0}) >> rshift) << Offset;
    mask = ~mask;
    current &= mask;

    mask = Bits == 64 ? ~word_t{0} : ((word_t{1} << Bits) - 1);
    current |= (static_cast<word_t>(value) & mask) << Offset;

    sim->state()[reg] = current;
  }

  PyMachine* sim;
};

RegisterFileRef PyMachine::getRegisters() {
  return RegisterFileRef(this);
}

#define REG64(r)                                                                                                       \
  def_property(                                                                                                        \
      #r, [](RegisterFileRef& r) { return r.getGPRegImpl<word_t, 64>(Register::r); },                                  \
      [](RegisterFileRef& r, word_t value) { r.setGPRegImpl<word_t, 64>(Register::r, value); })

#define REG32(name, r)                                                                                                 \
  def_property(                                                                                                        \
      #name, [](RegisterFileRef& r) { return r.getGPRegImpl<std::uint32_t, 32, 0>(Register::r); },                     \
      [](RegisterFileRef& r, std::uint32_t value) { r.setGPRegImpl<std::uint32_t, 32, 0, true>(Register::r, value); })

#define REG16(name, r)                                                                                                 \
  def_property(                                                                                                        \
      #name, [](RegisterFileRef& r) { return r.getGPRegImpl<std::uint16_t, 16>(Register::r); },                        \
      [](RegisterFileRef& r, std::uint16_t value) { r.setGPRegImpl<std::uint16_t, 16>(Register::r, value); })

#define REG8(name, r, offset)                                                                                          \
  def_property(                                                                                                        \
      #name, [](RegisterFileRef& r) { return r.getGPRegImpl<std::uint8_t, 8, offset>(Register::r); },                  \
      [](RegisterFileRef& r, std::uint8_t value) { r.setGPRegImpl<std::uint8_t, 8, offset>(Register::r, value); })

#define REG8L(name, r) REG8(name, r, 0)
#define REG8H(name, r) REG8(name, r, 8)

#define REGXMM(n)                                                                                                      \
  def_property_readonly("xmm" #n, [](RegisterFileRef& r) { return XMMRegister(r.sim, XmmRegister::xmm##n); })

#define EXCEPTION(name) py::register_exception<Raspsim##name>(m, #name, base_exc.ptr())

PYBIND11_MODULE(bindings, m) {
  m.doc() = "python binding for x86sim, a cycle-accurate x86 simulator based on PTLsim";

  py::class_<AddrRef>(m, "Address")
      .def("__add__", &AddrRef::operator+, "offset"_a, "Add an offset to the address")
      .def("__sub__", &AddrRef::operator-, "offset"_a, "Subtract an offset from the address")
      .def("__eq__", &AddrRef::operator==, "other"_a, "Check if two addresses are equal")
      .def("__ne__", &AddrRef::operator!=, "other"_a, "Check if two addresses are not equal")
      .def("__lt__", &AddrRef::operator<, "other"_a, "Check if the address is less than another address")
      .def("__le__", &AddrRef::operator<=, "other"_a, "Check if the address is less than or equal to another address")
      .def("__gt__", &AddrRef::operator>, "other"_a, "Check if the address is greater than another address")
      .def("__ge__", &AddrRef::operator>=, "other"_a,
           "Check if the address is greater than or equal to another address")
      .def("__hash__", [](const AddrRef& a) { return std::hash<address_t>{}(a.virtaddr); })
      .def("__int__", &AddrRef::operator address_t, "Get the address as an integer")
      .def("read", &AddrRef::read, "size"_a = 1, "Read data from the address")
      .def("write", &AddrRef::write, "value"_a, "Write data to the address");

  py::class_<MemImg>(m, "Memory")
      .def("__getitem__", &MemImg::getitem, "slice"_a, "Get a slice of memory")
      .def("__setitem__", &MemImg::setitem, "slice"_a, "data"_a, "Set a slice of memory");

  py::enum_<Prot>(m, "Prot")
      .value("READ", Prot::READ)
      .value("WRITE", Prot::WRITE)
      .value("EXEC", Prot::EXEC)
      .value("NONE", Prot::NONE)
      .value("RW", Prot::RW)
      .value("RX", Prot::RX)
      .value("RWX", Prot::RWX)
      .def("__or__", addProt, "Combine two protection flags")
      .def("__and__", hasProt, "Check if a protection flag is set");

  m.def("getProtFromELFSegment", &getProtFromELFSegment, "flags"_a,
        "Get the protection as Raspsim Prot from ELF segment flags");

  auto base_exc = py::register_exception<RaspsimException>(m, "RaspsimException");

  EXCEPTION(DivideException);
  EXCEPTION(DebugException);
  EXCEPTION(NMIException);
  EXCEPTION(BreakpointException);
  EXCEPTION(OverflowException);
  EXCEPTION(BoundsException);
  EXCEPTION(InvalidOpcodeException);
  EXCEPTION(FPUNotAvailException);
  EXCEPTION(DoubleFaultException);
  EXCEPTION(CoprocOverrunException);
  EXCEPTION(InvalidTSSException);
  EXCEPTION(SegNotPresentException);
  EXCEPTION(StackFaultException);
  EXCEPTION(GPFaultException);
  EXCEPTION(PageFaultException);
  EXCEPTION(SpuriousIntException);
  EXCEPTION(FPUException);
  EXCEPTION(UnalignedException);
  EXCEPTION(MachineCheckException);
  EXCEPTION(SSEException);

  py::class_<XMMRegister>(m, "XMMRegister", "A class to access the XMM registers of the virtual CPU")
      .def_property("sd", &XMMRegister::template getSingle<double>, &XMMRegister::template setSingle<double>)
      .def_property("ss", &XMMRegister::template getSingle<float>, &XMMRegister::template setSingle<float>)
      .def_property("pd", &XMMRegister::template getPacked<double, double>,
                    &XMMRegister::template setPacked<double, double>)
      .def_property("ps", &XMMRegister::template getPacked<float, float, float, float>,
                    &XMMRegister::template setPacked<float, float, float, float>)
      .def_property("chars",
                    &XMMRegister::template getPacked<char, char, char, char, char, char, char, char, char, char, char,
                                                     char, char, char, char, char>,
                    &XMMRegister::template setPacked<char, char, char, char, char, char, char, char, char, char, char,
                                                     char, char, char, char, char>);

  py::class_<RegisterFileRef>(m, "RegisterFile", "A class to access the registers of the virtual CPU")
      .def("__getitem__", &RegisterFileRef::getRegister, "regname"_a, "Get the value of a register")
      .def("__setitem__", &RegisterFileRef::setRegister, "regname"_a, "value"_a, "Set the value of a register")
      .REG64(rip)
      .REG64(rax)
      .REG64(rbx)
      .REG64(rcx)
      .REG64(rdx)
      .REG64(rsi)
      .REG64(rdi)
      .REG64(rbp)
      .REG64(rsp)
      .REG64(r8)
      .REG64(r9)
      .REG64(r10)
      .REG64(r11)
      .REG64(r12)
      .REG64(r13)
      .REG64(r14)
      .REG64(r15)
      .REG32(eax, rax)
      .REG32(ebx, rbx)
      .REG32(ecx, rcx)
      .REG32(edx, rdx)
      .REG32(esi, rsi)
      .REG32(edi, rdi)
      .REG32(ebp, rbp)
      .REG32(esp, rsp)
      .REG16(ax, rax)
      .REG16(bx, rbx)
      .REG16(cx, rcx)
      .REG16(dx, rdx)
      .REG16(si, rsi)
      .REG16(di, rdi)
      .REG16(bp, rbp)
      .REG16(sp, rsp)
      .REG8L(al, rax)
      .REG8L(bl, rbx)
      .REG8L(cl, rcx)
      .REG8L(dl, rdx)
      .REG8H(ah, rax)
      .REG8H(bh, rbx)
      .REG8H(ch, rcx)
      .REG8H(dh, rdx)
      .REGXMM(0)
      .REGXMM(1)
      .REGXMM(2)
      .REGXMM(3)
      .REGXMM(4)
      .REGXMM(5)
      .REGXMM(6)
      .REGXMM(7)
      .REGXMM(8)
      .REGXMM(9)
      .REGXMM(10)
      .REGXMM(11)
      .REGXMM(12)
      .REGXMM(13);

  py::class_<PyFs>(m, "Fs",
                   "Host-side handle to an in-memory guest filesystem (memfs).\n\n"
                   "Construct standalone to prepopulate a filesystem before creating a Machine "
                   "(pass it as Machine(memfs=fs)), share it between Machines, or obtain one from "
                   "Machine.fs. All paths are absolute or relative to the filesystem root. Errors "
                   "raise OSError subclasses (FileNotFoundError, FileExistsError, ...).")
      .def(py::init<>(), "Create a handle owning a fresh, empty filesystem")
      .def("listdir", &PyFs::listdir, "path"_a, "List the names in a directory (deterministic order)")
      .def("read_bytes", &PyFs::read_bytes, "path"_a, "Return the whole content of a file as bytes")
      .def("write_bytes", &PyFs::write_bytes, "path"_a, "data"_a,
           "Replace the content of a file (created if missing; parent directories must exist)")
      .def("mkdir", &PyFs::mkdir, "path"_a, "Create a directory (parent must exist)")
      .def("unlink", &PyFs::unlink, "path"_a, "Remove a file or symlink")
      .def("rmdir", &PyFs::rmdir, "path"_a, "Remove an empty directory")
      .def("rename", &PyFs::rename, "source"_a, "target"_a, "Rename/move with Linux rename(2) semantics")
      .def("symlink", &PyFs::symlink, "target"_a, "link_path"_a, "Create a symbolic link at link_path")
      .def("readlink", &PyFs::readlink, "path"_a, "Return the target of a symbolic link")
      .def("stat", &PyFs::stat, "path"_a, "follow_symlinks"_a = true,
           "Stat a path; returns a dict with st_mode/st_ino/st_nlink/st_size/st_atime/st_mtime/st_ctime. "
           "Timestamps are deterministic mutation ticks, not wall-clock time.")
      .def("truncate", &PyFs::truncate, "path"_a, "size"_a, "Truncate or zero-extend a file")
      .def("exists", &PyFs::exists, "path"_a, "Whether the path resolves to an existing node")
      .def("is_file", &PyFs::is_file, "path"_a, "Whether the path resolves to a regular file")
      .def("is_dir", &PyFs::is_dir, "path"_a, "Whether the path resolves to a directory")
      .def("is_symlink", &PyFs::is_symlink, "path"_a, "Whether the path itself is a symbolic link")
      // Two handles are equal when they view the same filesystem state.
      .def("__eq__",
           [](const PyFs& self, const py::object& other) {
             return py::isinstance<PyFs>(other) && self.state() == other.cast<const PyFs&>().state();
           })
      .def("__hash__", [](const PyFs& self) { return std::hash<const void*>{}(self.state().get()); });

  py::class_<PyMachine>(m, "Machine", "A class to interact with the simulator")
      .def(py::init<const char*, bool, bool, bool, bool, bool, py::object, py::object, py::object, py::object,
                    const char*, py::object>(),
           "logfile"_a = "", "sse"_a = true, "x87"_a = true, "perfect_cache"_a = false, "static_branchpred"_a = false,
           "glibc"_a = false, "stdin"_a = py::none(), "stdout"_a = py::none(), "stderr"_a = py::none(),
           "readlink"_a = py::none(), "core"_a = "ooo", "memfs"_a = py::none(),
           "Create a new Machine instance.\n\nSet glibc=True to enable the portable Linux "
           "glibc-startup syscalls: the malloc/free heap (brk + anonymous mmap/munmap/mremap) "
           "plus arch_prctl, set_tid_address, set_robust_list, rseq, prlimit64, uname, futex and "
           "the synthetic signal syscalls. Pass file-like objects for stdin/stdout/stderr to route "
           "guest read/write on fds 0/1/2 to Python (e.g. io.BytesIO); host fds are never used. Pass "
           "readlink=callable to resolve guest readlink()/readlinkat() calls (path: str -> str | "
           "None), e.g. for /proc/self/exe; without it readlink returns -ENOSYS.\n\n"
           "Set memfs=True (or pass an x86sim Fs object to share/prepopulate one) to give the guest "
           "an in-memory Linux-like filesystem: open/read/write/stat/getdents64/... work against a "
           "portable sandbox that the host can inspect through Machine.fs. No fd number is special: "
           "stdio stays with the stdin/stdout/stderr streams unless the guest dup2()s over them or "
           "map_fd() binds a memfs file there. memfs is mutually exclusive with readlink.\n\n"
           "core selects the CPU model: \"ooo\" (default, out-of-order) or \"seq\" (sequential). "
           "The sequential core is slower but handles unaligned memory accesses correctly; the "
           "out-of-order core currently cannot, so running glibc (which uses unaligned SSE in its "
           "string routines) requires core=\"seq\".")
      .def_property_readonly("registers", &PyMachine::getRegisters, "Get the register file")
      .def("run", &PyMachine::run, "Run the simulator for a number of instructions",
           "ninstr"_a = static_cast<unsigned long long>(-1))
      .def_property_readonly("cycles", &PyMachine::cycles, "Get the number of cycles")
      .def_property_readonly("instructions", &PyMachine::instructions, "Get the number of instructions")
      .def("memmap", &PyMachine::memmap, "start"_a, "prot"_a, "length"_a = 0, "data"_a = py::bytes(),
           "Map a range of memory to the virtual address space of the "
           "simulator.\n\nMaps data from `data` into memory and fills the "
           "rest with zeros if `length` is greater than the size of `data`. "
           "If `length` is 0, the size of `data` will be used as length.")
      .def_property_readonly(
          "memimg", [](PyMachine& sim) { return MemImg(sim); }, "Get a memory image object")
      .def_property_readonly("fs", &PyMachine::fs,
                             "Host-side Fs handle to the guest's in-memory filesystem (memfs=True only)")
      .def("map_fd", &PyMachine::map_fd, "fd"_a, "path"_a, "mode"_a = "r",
           "Bind a memfs file to an arbitrary guest fd number (including 0/1/2, replacing a "
           "configured stream). mode follows open(): \"r\", \"w\", \"a\", optionally with \"+\".")
      .def("__str__", &PyMachine::str, "Get the string representation of the current state of the simulator");
}
