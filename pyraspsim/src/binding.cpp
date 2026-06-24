#include <pybind11/pybind11.h>

#include "x86sim-support/cpuid.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <format>
#include <mutex>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
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

class PyRaspsim;

class AddrRef {
public:
  AddrRef(PyRaspsim* sim, address_t virtaddr) : virtaddr(virtaddr), sim(sim) {}
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
  PyRaspsim* sim;
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
  class name : public RaspsimException {                                                                              \
  public:                                                                                                              \
    name(const RaspsimException& e) : RaspsimException(e) {}                                                          \
    name(RaspsimException&& e) : RaspsimException(std::move(e)) {}                                                    \
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

// raspsim is a register/memory-level simulator: an int 0x80 signals the guest is
// done; any other syscall is unsupported and surfaces as an exception.
class PyHost : public x86sim::HostCallbacks {
public:
  x86sim::SyscallResult syscall(x86sim::Machine&, RegisterFile&, x86sim::SyscallKind kind) override {
    using x86sim::StopReason;
    if (kind == x86sim::SyscallKind::int80)
      return {.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};
    return {.reason = StopReason::unsupported_syscall, .continue_execution = false,
            .message = "Syscall not supported"};
  }

  x86sim::CpuidResult cpuid(x86sim::Machine&, RegisterFile&, x86sim::CpuidRequest request) noexcept override {
    return x86sim::defaults::default_cpuid(request);
  }
};

class PyRaspsim {
public:
  PyRaspsim(const char* logfile, bool sse, bool x87, bool perfect_cache, bool static_branchpred) {
    if (!lock.try_lock())
      throw py::value_error("Only one instance of Raspsim can be used at a time");

    x86sim::Options options;
    options.sse = sse;
    options.x87 = x87;
    options.debug.perfect_cache = perfect_cache;
    options.debug.static_branchpred = static_branchpred;
    options.log.log_filename = logfile;

    machine = std::make_unique<x86sim::Machine>(host, options);
  }

  ~PyRaspsim() { lock.unlock(); }

  static constexpr word_t getPageSize() { return x86sim::Machine::kPageSize; }

  x86sim::Machine& m() { return *machine; }

  RegisterFileRef getRegisters();

  AddrRef memmap(address_t start, Prot prot, word_t length, py::bytes data);

  std::size_t cycles() { return machine->stats().cycles; }
  std::size_t instructions() { return machine->stats().instructions; }

  std::string str() { return std::format("{}", machine->register_file(0)); }

  void run(unsigned long long ninstr);

  PyHost host;
  std::unique_ptr<x86sim::Machine> machine;
  static std::mutex lock;
};

std::mutex PyRaspsim::lock;

void AddrRef::write(py::bytes&& bts) {
  std::string mem{std::move(bts)};
  auto bytes = std::as_bytes(std::span(mem.data(), mem.size()));
  if (auto r = sim->m().write_memory(virtaddr, bytes); !r)
    throw py::value_error(std::format("Trying to write to unmapped memory at {:x}: {}", virtaddr, r.error()));
}

py::bytes AddrRef::read(word_t size) {
  std::string out(size, '\0');
  auto bytes = std::as_writable_bytes(std::span(out.data(), out.size()));
  if (auto r = sim->m().read_memory_into(virtaddr, bytes); !r)
    throw py::value_error(std::format("Trying to read from unmapped memory at {:x}: {}", virtaddr, r.error()));
  return {std::move(out)};
}

AddrRef PyRaspsim::memmap(address_t start, Prot prot, word_t length, py::bytes data) {
  std::string bytes = std::move(data);
  if (bytes.size() > 0) {
    if (length > 0 && bytes.size() > length)
      throw py::value_error("Data size must be less than or equal to length");
    length = std::max(static_cast<word_t>(length), static_cast<word_t>(bytes.size()));
  }

  if (length == 0)
    throw py::value_error("Cannot map zero length data");

  word_t offset = start % getPageSize();
  if (auto r = machine->map(start - offset, length + offset, toProtection(prot)); !r)
    throw py::value_error(std::format("Cannot map memory at {:x}: {}", start, r.error()));

  if (!bytes.empty()) {
    auto raw = std::as_bytes(std::span(bytes.data(), bytes.size()));
    if (auto r = machine->write_memory(start, raw); !r)
      throw py::value_error(std::format("Cannot write mapped memory at {:x}: {}", start, r.error()));
  }
  return AddrRef(this, start);
}

#define THROW(i, cls)                                                                                                  \
  case i:                                                                                                              \
    throw cls(std::move(exc));

void PyRaspsim::run(unsigned long long ninstr) {
  x86sim::RunOptions run_options;
  if (ninstr != static_cast<unsigned long long>(-1))
    run_options.instruction_limit = ninstr;

  x86sim::RunResult result = machine->run(run_options);

  using x86sim::StopReason;
  switch (result.reason) {
  case StopReason::guest_exit:
    return;
  case StopReason::instruction_limit:
    throw py::stop_iteration("Reached instruction limit");
  case StopReason::x86_exception: {
    RaspsimException exc = result.x86_exception ? RaspsimException(*result.x86_exception)
                                                : RaspsimException("Unknown x86 exception");
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
      {"rax", rax}, {"rcx", rcx}, {"rdx", rdx}, {"rbx", rbx}, {"rsp", rsp}, {"rbp", rbp}, {"rsi", rsi}, {"rdi", rdi},
      {"r8", r8},   {"r9", r9},   {"r10", r10}, {"r11", r11}, {"r12", r12}, {"r13", r13}, {"r14", r14}, {"r15", r15},
      {"rip", rip}, {"flags", flags},
  };
  for (auto [reg_name, reg] : names) {
    if (name == reg_name)
      return reg;
  }
  return std::nullopt;
}

class XMMRegister {
public:
  XMMRegister(PyRaspsim* sim, XmmRegister reg) : sim(sim), reg(reg) {}

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
    XmmValue value = sim->m().register_file(0)[reg];
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
    sim->m().register_file(0)[reg] = XmmValue{raw[0], raw[1]};
  }

  template<typename T>
  void setSingle(T value) {
    setPacked<T>(std::tuple<T>{value});
  }

private:
  PyRaspsim* sim;
  XmmRegister reg;
};

class MemImg {
public:
  MemImg(PyRaspsim& sim) : sim(&sim) {}

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

  PyRaspsim* sim;
};

class RegisterFileRef {
  friend class XMMRegister;

public:
  RegisterFileRef() = delete;
  RegisterFileRef(PyRaspsim* sim) : sim(sim) {}

  word_t getRegister(std::string&& regname) {
    auto reg = registerFromName(regname);
    if (!reg)
      throw py::value_error(std::string("Invalid register name '") + regname + "'");
    return sim->m().register_file(0)[*reg];
  }

  void setRegister(std::string&& regname, word_t value) {
    auto reg = registerFromName(regname);
    if (!reg)
      throw py::value_error(std::string("Invalid register name '") + regname + "'");
    sim->m().register_file(0)[*reg] = value;
  }

  template<typename T, int Bits, int Offset = 0>
  inline T getGPRegImpl(Register reg) {
    static_assert(Bits <= 64, "Register size must be less than or equal to 64 bits");
    static_assert(Offset % 8 == 0, "Offset must be a multiple of 8");
    static_assert(Bits % 8 == 0, "Bits must be a multiple of 8");

    word_t mask = Bits == 64 ? ~word_t{0} : ((word_t{1} << Bits) - 1);
    return (sim->m().register_file(0)[reg] >> Offset) & mask;
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

    word_t current = sim->m().register_file(0)[reg];
    int rshift = ZeroHigh ? 0 : 64 - Bits;

    // mask contains 1s in the bits that are not part of the register that
    // should be modified
    word_t mask = ((~word_t{0}) >> rshift) << Offset;
    mask = ~mask;
    current &= mask;

    mask = Bits == 64 ? ~word_t{0} : ((word_t{1} << Bits) - 1);
    current |= (static_cast<word_t>(value) & mask) << Offset;

    sim->m().register_file(0)[reg] = current;
  }

  PyRaspsim* sim;
};

RegisterFileRef PyRaspsim::getRegisters() {
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

PYBIND11_MODULE(core, m) {
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

  py::class_<PyRaspsim>(m, "Core", "A class to interact with the simulator")
      .def(py::init<const char*, bool, bool, bool, bool>(), "logfile"_a = "/dev/null", "sse"_a = true, "x87"_a = true,
           "perfect_cache"_a = false, "static_branchpred"_a = false, "Create a new Raspsim instance")
      .def_property_readonly("registers", &PyRaspsim::getRegisters, "Get the register file")
      .def("run", &PyRaspsim::run, "Run the simulator for a number of instructions",
           "ninstr"_a = static_cast<unsigned long long>(-1))
      .def_property_readonly("cycles", &PyRaspsim::cycles, "Get the number of cycles")
      .def_property_readonly("instructions", &PyRaspsim::instructions, "Get the number of instructions")
      .def("memmap", &PyRaspsim::memmap, "start"_a, "prot"_a, "length"_a = 0, "data"_a = py::bytes(),
           "Map a range of memory to the virtual address space of the "
           "simulator.\n\nMaps data from `data` into memory and fills the "
           "rest with zeros if `length` is greater than the size of `data`. "
           "If `length` is 0, the size of `data` will be used as length.")
      .def_property_readonly(
          "memimg", [](PyRaspsim& sim) { return MemImg(sim); }, "Get a memory image object")
      .def("__str__", &PyRaspsim::str, "Get the string representation of the current state of the simulator");
}
