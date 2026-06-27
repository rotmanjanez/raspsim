#include <pybind11/pybind11.h>

#include "x86sim-support/cpuid.hpp"
#include "x86sim-support/syscall-linux.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <format>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>

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

private:
  // Cap a single read/write transfer so a bogus guest count cannot trigger a
  // huge allocation; oversized requests become a (legal) short transfer. This
  // mirrors the spirit of the bounds checks in the syscall support library.
  static constexpr word_t kMaxTransfer = word_t{1} << 31; // 2 GiB

  // x86-64 syscall numbers handled directly here (routed to Python), and the
  // negated Linux errno values returned to the guest on the failure paths.
  static constexpr word_t kSyscallReadlink = 89;
  static constexpr word_t kSyscallReadlinkat = 267;
  static constexpr std::int64_t kLinuxEnoent = 2;
  static constexpr std::int64_t kLinuxEnosys = 38;

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
           x86sim::linux_syscalls::SysMremap{} | x86sim::linux_syscalls::SysArchPrctl{} |
           x86sim::linux_syscalls::SysSetTidAddress{} | x86sim::linux_syscalls::SysSetRobustList{} |
           x86sim::linux_syscalls::SysRseq{} | x86sim::linux_syscalls::SysPrlimit64{} |
           x86sim::linux_syscalls::SysUname{} | x86sim::linux_syscalls::SysGetIdentity{} |
           x86sim::linux_syscalls::SysFutex{} | x86sim::linux_syscalls::SysSignals{} |
           x86sim::linux_syscalls::SysGetrandom{}) glibc_syscalls =
      x86sim::linux_syscalls::SysBrk{} | x86sim::linux_syscalls::SysMmap{} | x86sim::linux_syscalls::SysMunmap{} |
      x86sim::linux_syscalls::SysMremap{} | x86sim::linux_syscalls::SysArchPrctl{} |
      x86sim::linux_syscalls::SysSetTidAddress{} | x86sim::linux_syscalls::SysSetRobustList{} |
      x86sim::linux_syscalls::SysRseq{} | x86sim::linux_syscalls::SysPrlimit64{} | x86sim::linux_syscalls::SysUname{} |
      x86sim::linux_syscalls::SysGetIdentity{} | x86sim::linux_syscalls::SysFutex{} |
      x86sim::linux_syscalls::SysSignals{} | x86sim::linux_syscalls::SysGetrandom{};
};

class PyMachine {
public:
  PyMachine(const char* logfile, bool sse, bool x87, bool perfect_cache, bool static_branchpred, bool glibc,
            py::object stdin_obj, py::object stdout_obj, py::object stdout_err, py::object readlink_cb,
            const char* core) {
    host.enable_glibc = glibc;
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

    machine = std::make_unique<x86sim::Machine>(host, options);
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

  PyHost host;
  x86sim::CpuState cpu_state;
  x86sim::AddressSpace address_space;
  std::unique_ptr<x86sim::Machine> machine;

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

  py::class_<PyMachine>(m, "Machine", "A class to interact with the simulator")
      .def(py::init<const char*, bool, bool, bool, bool, bool, py::object, py::object, py::object, py::object,
                    const char*>(),
           "logfile"_a = "/dev/null", "sse"_a = true, "x87"_a = true, "perfect_cache"_a = false,
           "static_branchpred"_a = false, "glibc"_a = false, "stdin"_a = py::none(), "stdout"_a = py::none(),
           "stderr"_a = py::none(), "readlink"_a = py::none(), "core"_a = "ooo",
           "Create a new Machine instance.\n\nSet glibc=True to enable the portable Linux "
           "glibc-startup syscalls: the malloc/free heap (brk + anonymous mmap/munmap/mremap) "
           "plus arch_prctl, set_tid_address, set_robust_list, rseq, prlimit64, uname, futex and "
           "the synthetic signal syscalls. Pass file-like objects for stdin/stdout/stderr to route "
           "guest read/write on fds 0/1/2 to Python (e.g. io.BytesIO); host fds are never used. Pass "
           "readlink=callable to resolve guest readlink()/readlinkat() calls (path: str -> str | "
           "None), e.g. for /proc/self/exe; without it readlink returns -ENOSYS.\n\n"
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
      .def("__str__", &PyMachine::str, "Get the string representation of the current state of the simulator");
}
