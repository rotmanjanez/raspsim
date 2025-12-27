#include <elf.h>
#include <pybind11/pybind11.h>
#include <raspsim-hwsetup.h>
#include <setjmp.h>
#include <sys/mman.h>

#include <algorithm>
#include <mutex>
#include <registers.def>
#include <sstream>
#include <string>

namespace py = pybind11;
using namespace py::literals;

enum class Prot {
    READ = PROT_READ,
    WRITE = PROT_WRITE,
    EXEC = PROT_EXEC,
    NONE = PROT_NONE,
    RW = PROT_READ | PROT_WRITE,
    RX = PROT_READ | PROT_EXEC,
    RWX = PROT_READ | PROT_WRITE | PROT_EXEC
};

std::ostream &operator<<(std::ostream &os, Prot p) {
    switch (p) {
        case Prot::READ:
            return os << "READ";
        case Prot::WRITE:
            return os << "WRITE";
        case Prot::EXEC:
            return os << "EXEC";
        case Prot::NONE:
            return os << "NONE";
        case Prot::RW:
            return os << "RW";
        case Prot::RX:
            return os << "RX";
        case Prot::RWX:
            return os << "RWX";
    }

    return os;
}

bool hasProt(Prot p, Prot q) {
    return (static_cast<int>(p) & static_cast<int>(q)) == static_cast<int>(q);
}

Prot addProt(Prot p, Prot q) {
    return static_cast<Prot>(static_cast<int>(p) | static_cast<int>(q));
}

Prot getProtFromELFSegment(int flags) {
    int prot = 0;

    if (flags & PF_R) {
        prot |= static_cast<int>(Prot::READ);
    }
    if (flags & PF_W) {
        prot |= static_cast<int>(Prot::WRITE);
    }

    if (flags & PF_X) {
        prot |= static_cast<int>(Prot::EXEC);
    }

    return static_cast<Prot>(prot);
}

class AddrRef {
   public:
    AddrRef(Raspsim *sim, Waddr virtaddr) : virtaddr(virtaddr), sim(sim) {}
    AddrRef() : virtaddr(0), sim(nullptr) {}

    operator Waddr() const { return virtaddr; }

    AddrRef operator+(Waddr offset) { return AddrRef{sim, virtaddr + offset}; }
    AddrRef operator-(Waddr offset) { return AddrRef{sim, virtaddr - offset}; }

    bool operator==(const AddrRef &other) const {
        return addr() == other.addr();
    }
    bool operator!=(const AddrRef &other) const {
        return addr() != other.addr();
    }
    bool operator<(const AddrRef &other) const { return addr() < other.addr(); }
    bool operator<=(const AddrRef &other) const {
        return addr() <= other.addr();
    }
    bool operator>(const AddrRef &other) const { return addr() > other.addr(); }
    bool operator>=(const AddrRef &other) const {
        return addr() >= other.addr();
    }

    void write(py::bytes &&bts) {
        std::string mem{std::move(bts)};

        W64 ps = Raspsim::getPageSize();
        W64 nprocessed = 0;

        while (nprocessed != mem.size()) {
            W64 chunksize = std::min(
                {mem.size() - nprocessed, (W64)ps, ps - (virtaddr % (ps - 1))});
            void *addr = sim->page_virt_to_mapped(virtaddr + nprocessed);
            if (!addr) {
                throw py::value_error("Trying to write to unmapped memory");
            }
            memcpy(addr, mem.data() + nprocessed, chunksize);
            nprocessed += chunksize;
        }
    }

    py::bytes read(W64 size = 1) {
        std::string b = "";
        size_t ps = Raspsim::getPageSize();
        W64 nprocessed = 0;

        while (size) {
            Waddr effvirtaddr = virtaddr + nprocessed;
            W64 chunksize = std::min({(W64)size, (W64)ps, ps - (effvirtaddr & (ps - 1))});
            void *addr = sim->page_virt_to_mapped(effvirtaddr);
            if (!addr) {
                std::stringstream ss;
                ss << "Trying to read from unmapped memory at " << std::hex << effvirtaddr << std::dec;
                throw py::value_error(strdup(ss.str().c_str()));
            }
            b += std::string((char *)addr, chunksize);
            nprocessed += chunksize;
            size -= chunksize;
        }

        return {std::move(b)};
    }

    void *addr() const { return sim->page_virt_to_mapped(virtaddr); }

    Waddr virtaddr;

   protected:
    Raspsim *sim;
};

extern bool ensureMachineInitialized(PTLsimMachine &m, const char *machinename);
extern void simulateInitializedMachine(PTLsimMachine &m);

class RaspsimException : public std::exception {
   public:
    RaspsimException(byte exception, W32 errorcode, Waddr virtaddr)
        : exception(exception),
          errorcode(errorcode),
          virtaddr(virtaddr),
          msg("") {
        char *m = Raspsim::formatException(exception, errorcode, virtaddr);
        msg += m;
        free(m);
    }
    RaspsimException(byte exception, W32 errorcode, Waddr virtaddr,
                     const Context &ctx)
        : RaspsimException(exception, errorcode, virtaddr) {
        msg += "\n";
        msg += Raspsim::formatContext(ctx);
    }

    RaspsimException(const char *msg)
        : exception(-1), errorcode(0), virtaddr(0), msg(strdup(msg)) {}

    byte getException() const { return exception; }
    W32 getErrorCode() const { return errorcode; }
    Waddr getVirtAddr() const { return virtaddr; }

    const char *what() const noexcept override { return msg.c_str(); }

    std::pair<std::string, unsigned> loc{};

   private:
    byte exception;
    W32 errorcode;
    Waddr virtaddr;
    std::string msg;
};

#define RASPSIM_INHERIT_EXCEPTION(name)                                \
    class name : public RaspsimException {                             \
       public:                                                         \
        name(const RaspsimException &e) : RaspsimException(e) {}       \
        name(RaspsimException &&e) : RaspsimException(std::move(e)) {} \
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

class RegisterFile;

class PyRaspsim : public Raspsim {
   public:
    PyRaspsim(const char *logfile) : Raspsim() {
        if (lock.try_lock()) {
            setLogfile(logfile);
        } else {
            throw py::value_error(
                "Only one instance of Raspsim can be used at a time");
        }
    }

    ~PyRaspsim() { lock.unlock(); }

    RegisterFile getRegisters();

    AddrRef memmap(Waddr start, Prot prot, Waddr length, py::bytes data);

    std::size_t cycles() { return Raspsim::cycles(); }
    std::size_t instructions() { return Raspsim::instructions(); }

    std::string str() {
        const char *c = Raspsim::formatContext(Raspsim::getContext());
        std::string s(c);
        free((void *)c);
        return s;
    }

    void run(unsigned long long ninstr = (unsigned long long)-1);

    static jmp_buf simexit;
    static RaspsimException X86Exception;
    static std::mutex lock;
};

std::mutex PyRaspsim::lock;

AddrRef PyRaspsim::memmap(Waddr start, Prot prot, W64 length, py::bytes data) {
    std::string bytes = std::move(data);
    byte *raw = nullptr;
    if (bytes.size() > 0) {
        raw = (byte *)bytes.data();
        if (length > 0 && bytes.size() > length) {
            throw py::value_error(
                "Data size must be less than or equal to length");
        }
        length = std::max((std::size_t)length, bytes.size());
    }

    if (length == 0) {
        throw py::value_error("Cannot map zero length data");
    }

    W64 offset = start % getPageSize();

    map(start, length + offset, static_cast<int>(prot));

    if (raw) {
        const W64 pg = getPageSize();
        W64 nprocessed = 0;

        while (nprocessed < bytes.size()) {
            Waddr virtaddr = start + nprocessed;
            W64 offset = virtaddr & (pg - 1);
            W64 chunksize = std::min(bytes.size() - nprocessed, pg - offset);
            
            void *dst = page_virt_to_mapped(virtaddr);

            memcpy(dst, raw + nprocessed, chunksize);
            nprocessed += chunksize;
        }
    }
    return AddrRef(this, start);
}

#define THROW(i, cls)                       \
    case i:                                 \
        throw cls(std::move(X86Exception)); \
        break;

void PyRaspsim::run(unsigned long long ninstr) {
    setTimeout(ninstr);
    PTLsimMachine *machine = getMachine();

    if (!machine) {
        throw py::value_error(std::string("Cannot find core named '") +
                              getCoreName() + "'");
    }

    if (ensureMachineInitialized(*machine, getCoreName())) {
        throw std::runtime_error(
            std::string("Cannot initialize core model for '") + getCoreName() +
            "'");
    }

    if (setjmp(simexit)) {
        switch (X86Exception.getException()) {
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
                throw X86Exception;
        }
    } else {
        simulateInitializedMachine(*machine);
        if (instructions() >= ninstr) {
            throw py::stop_iteration("Reached instruction limit");
        }
    }
}

jmp_buf PyRaspsim::simexit;
RaspsimException PyRaspsim::X86Exception = {0, 0, 0};

#define XMM(number)                                                            \
    void set_xmml##number(W64 value) {                                         \
        sim->setRegisterValue(REG_xmml##number, value);                        \
    }                                                                          \
    void set_xmmh##number(W64 value) {                                         \
        sim->setRegisterValue(REG_xmmh##number, value);                        \
    }                                                                          \
    W64 get_xmml##number() { return sim->getRegisterValue(REG_xmml##number); } \
    W64 get_xmmh##number() { return sim->getRegisterValue(REG_xmmh##number); }

class XMMRegister {
   public:
    XMMRegister(PyRaspsim *sim, int reg) : sim(sim), reg(reg) {
        assert(!(reg & 1) &&
               "XMM Register must always be indexed via the lower 64 Bits");
    }

    template <typename T, typename... Ts>
    std::tuple<T, Ts...> getPacked() {
        constexpr int NItems = sizeof...(Ts) + 1;
        static_assert(sizeof(T) <= 8,
                      "XMM Register can only be cast to vector types with an "
                      "element type size less than or equal to 64 bits");
        static_assert(
            (NItems * sizeof(T) == 16) || (NItems == 1 && sizeof(T) <= 8),
            "XMM Register can only be cast to vector types with size equal to "
            "128 bits or scalar types with a size less than or equal to 64 "
            "bits");

        std::tuple<T, Ts...> result;
        std::array<W64, 2> value;

        value[0] = sim->getRegisterValue(reg);
        value[1] = sim->getRegisterValue(reg + 1);

        memccpy(static_cast<void *>(&std::get<0>(result)), value.data(), NItems,
                sizeof(T));
        return result;
    }

    template <typename T>
    T getSingle() {
        return std::get<0>(getPacked<T>());
    }

    template <typename T, typename... Ts>
    void setPacked(std::tuple<T, Ts...> value) {
        constexpr int NItems = sizeof...(Ts) + 1;
        static_assert(sizeof(T) <= 8,
                      "XMM Register can only be cast to vector types with an "
                      "element type size less than or equal to 64 bits");
        static_assert(
            (NItems * sizeof(T) == 16) || (NItems == 1 && sizeof(T) <= 8),
            "XMM Register can only be cast to vector types with size equal to "
            "128 bits or scalar types with a size less than or equal to 64 "
            "bits");

        std::array<W64, 2> result = {0, 0};
        memccpy(result.data(), static_cast<void *>(&std::get<0>(value)), NItems,
                sizeof(T));

        sim->setRegisterValue(reg, result[0]);
        sim->setRegisterValue(reg + 1, result[1]);
    }

    template <typename T>
    void setSingle(T value) {
        setPacked<T>({value});
    }

   private:
    PyRaspsim *sim;
    int reg;
};

class MemImg {
   public:
    MemImg(PyRaspsim &sim) : sim(&sim) {}
    py::bytes getitem(const py::slice& s) const {
        auto startstop = validateSlice(s);
        
        return AddrRef(sim, startstop.first).read(startstop.second - startstop.first);
    }

    void setitem(const py::slice& s, py::bytes data) {
        auto startstop = validateSlice(s);
        std::string raw{std::move(data)};

        if (raw.size() != startstop.second - startstop.first) {
            throw py::value_error("Data size must match slice size");
        }

        AddrRef(sim, startstop.first).write(std::move(data));
    }

   private:

   std::pair<Waddr, W64> validateSlice(const py::slice& s) const {
        if (!py::isinstance<py::none>(s.attr("step"))) {
            throw py::value_error("Step is not supported");
        }

        return {s.attr("start").cast<Waddr>(), s.attr("stop").cast<W64>()};
    }

    PyRaspsim *sim;
};

class RegisterFile {
    friend class XMMRegister;

   public:
    RegisterFile() = delete;

    RegisterFile(PyRaspsim *sim) : sim(sim) {}

    RegisterFile(const RegisterFile &) = default;
    RegisterFile(RegisterFile &&) = default;
    RegisterFile &operator=(const RegisterFile &) = delete;
    RegisterFile &operator=(RegisterFile &&) = delete;

    W64 getRegister(std::string &&regname) {
        return sim->getRegisterValue(getRegisterIdx(std::move(regname)));
    }

    void setRegister(std::string &&regname, W64 value) {
        sim->setRegisterValue(getRegisterIdx(std::move(regname)), value);
    }

    template <typename T, int Bits, int Offset = 0>
    inline T getGPRegImpl(int reg) {
        static_assert(Bits <= 64,
                      "Register size must be less than or equal to 64 bits");
        static_assert(Offset % 8 == 0, "Offset must be a multiple of 8");
        static_assert(Bits % 8 == 0, "Bits must be a multiple of 8");

        W64 mask = Bits == 64 ? -1 : (((W64)1 << Bits) - 1);
        return (sim->getRegisterValue(reg) >> Offset) & mask;
    }

    template <typename T, unsigned Bits, unsigned Offset = 0,
              bool ZeroHigh = false>
    inline void setGPRegImpl(int reg, T value) {
        static_assert(Bits <= 64,
                      "Register size must be less than or equal to 64 bits");
        static_assert(Offset + Bits <= 64,
                      "Offset + Bits must be less than or equal to 64");
        static_assert(Offset % 8 == 0, "Offset must be a multiple of 8");
        static_assert(Bits % 8 == 0, "Bits must be a multiple of 8");
        static_assert(
            !(ZeroHigh && Bits != 32 && Offset != 0),
            "ZeroHigh can only be used with 32-bit registers at offset 0");
        static_assert(!(Offset > 0 && Bits != 8),
                      "Offset can only be used with 8-bit registers");

        W64 current = sim->getRegisterValue(reg);
        int rshift = ZeroHigh ? 0 : 64 - Bits;

        // mask containts 1s in the bits that are not part of the register that
        // should be modified
        W64 mask = (((W64)-1) >> rshift) << Offset;

        mask = ~mask;
        current &= mask;

        mask = Bits == 64 ? -1 : (((W64)1 << Bits) - 1);
        current |= (value & mask) << Offset;

        sim->setRegisterValue(reg, current);
    }

   private:
    unsigned getRegisterIdx(std::string &&regname) {
        int idx = sim->getRegisterIndex(regname.c_str());
        if (idx < 0) {
            throw py::value_error(std::string("Invalid register name '") +
                                  regname + "'");
        }
        return idx;
    }

   public:
    PyRaspsim *sim;
};

RegisterFile PyRaspsim::getRegisters() { return RegisterFile(this); }

#define REG64(r)                                                              \
    def_property(                                                             \
        #r, [](RegisterFile &r) { return r.getGPRegImpl<W64, 64>(REG_##r); }, \
        [](RegisterFile &r, W64 value) {                                      \
            r.setGPRegImpl<W64, 64>(REG_##r, value);                          \
        })

#define REG32(name, r)                                                       \
    def_property(                                                            \
        #name,                                                               \
        [](RegisterFile &r) { return r.getGPRegImpl<W32, 32, 0>(REG_##r); }, \
        [](RegisterFile &r, W32 value) {                                     \
            r.setGPRegImpl<W32, 32, 0, true>(REG_##r, value);                \
        })

#define REG16(name, r)                                                    \
    def_property(                                                         \
        #name,                                                            \
        [](RegisterFile &r) { return r.getGPRegImpl<W16, 16>(REG_##r); }, \
        [](RegisterFile &r, W16 value) {                                  \
            r.setGPRegImpl<W16, 16>(REG_##r, value);                      \
        })

#define REG8(name, r, offset)                              \
    def_property(                                          \
        #name,                                             \
        [](RegisterFile &r) {                              \
            return r.getGPRegImpl<W8, 8, offset>(REG_##r); \
        },                                                 \
        [](RegisterFile &r, W8 value) {                    \
            r.setGPRegImpl<W8, 8, offset>(REG_##r, value); \
        })

#define REG8L(name, r) REG8(name, r, 0)

#define REG8H(name, r) REG8(name, r, 8)

#define REGXMM(n)                                         \
    def_property_readonly("xmm" #n, [](RegisterFile &r) { \
        return XMMRegister(r.sim, REG_xmml##n);           \
    })

#define EXCEPTION(name) py::register_exception<Raspsim##name>(m, #name)

PYBIND11_MODULE(core, m) {
    m.doc() =
        "python binding for raspsim, a cycle-accurate x86 simulator based on "
        "PTLsim";

    py::class_<AddrRef>(m, "Address")
        .def("__add__", &AddrRef::operator+, "offset"_a,
             "Add an offset to the address")
        .def("__sub__", &AddrRef::operator-, "offset"_a,
             "Subtract an offset from the address")
        .def("__eq__", &AddrRef::operator==, "other"_a,
             "Check if two addresses are equal")
        .def("__ne__", &AddrRef::operator!=, "other"_a,
             "Check if two addresses are not equal")
        .def("__lt__", &AddrRef::operator<, "other"_a,
             "Check if the address is less than another address")
        .def("__le__", &AddrRef::operator<=, "other"_a,
             "Check if the address is less than or equal to another address")
        .def("__gt__", &AddrRef::operator>, "other"_a,
             "Check if the address is greater than another address")
        .def("__ge__", &AddrRef::operator>=, "other"_a,
             "Check if the address is greater than or equal to another address")
        .def("__hash__",
             [](const AddrRef &a) { return std::hash<Waddr>{}(a.virtaddr); })
        .def("__int__", &AddrRef::operator Waddr, "Get the address as an integer")
        .def("read", &AddrRef::read, "size"_a = 1, "Read data from the address")
        .def("write", &AddrRef::write, "value"_a, "Write data to the address");

    py::class_<MemImg>(m, "Memory")
        .def("__getitem__", &MemImg::getitem, "slice"_a,
             "Get a slice of memory")
        .def("__setitem__", &MemImg::setitem, "slice"_a, "data"_a,
                "Set a slice of memory");

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

    py::register_exception<RaspsimException>(m, "RaspsimException");

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

    py::class_<XMMRegister>(
        m, "XMMRegister",
        "A class to access the XMM registers of the virtual CPU")
        .def_property("sd", &XMMRegister::template getSingle<double>,
                      &XMMRegister::template setSingle<double>)
        .def_property("ss", &XMMRegister::template getSingle<float>,
                      &XMMRegister::template getSingle<float>)
        .def_property("pd", &XMMRegister::template getPacked<double, double>,
                      &XMMRegister::template setPacked<double, double>)
        .def_property(
            "ps", &XMMRegister::template getPacked<float, float, float, float>,
            &XMMRegister::template setPacked<float, float, float, float>)
        .def_property(
            "chars",
            &XMMRegister::template getPacked<char, char, char, char, char, char,
                                             char, char, char, char, char, char,
                                             char, char, char, char>,
            &XMMRegister::template setPacked<char, char, char, char, char, char,
                                             char, char, char, char, char, char,
                                             char, char, char, char>);

    py::class_<RegisterFile>(
        m, "RegisterFile", "A class to access the registers of the virtual CPU")
        .def("__getitem__", &RegisterFile::getRegister, "regname"_a,
             "Get the value of a register")
        .def("__setitem__", &RegisterFile::setRegister, "regname"_a, "value"_a,
             "Set the value of a register")
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
        .def(py::init<const char *>(), "logfile"_a = "/dev/null",
             "Create a new Raspsim instance")
        .def_property_readonly("registers", &PyRaspsim::getRegisters,
                               "Get the register file")
        .def("run", &PyRaspsim::run,
             "Run the simulator for a number of instructions",
             "ninstr"_a = (unsigned long long)-1)
        .def("disableSSE", &PyRaspsim::disableSSE, "Disable SSE")
        .def("disableX87", &PyRaspsim::disableX87, "Disable X87")
        .def("enablePerfectCache", &PyRaspsim::enablePerfectCache,
             "Enable perfect cache")
        .def("enableStaticBranchPrediction",
             &PyRaspsim::enableStaticBranchPrediction,
             "Enable static branch prediction")
        .def_property_readonly("cycles", &PyRaspsim::cycles,
                               "Get the number of cycles")
        .def_property_readonly("instructions", &PyRaspsim::instructions,
                               "Get the number of instructions")
        .def("memmap", &PyRaspsim::memmap,
             "Map a range of memory to the virtual address space of the "
             "simulator",
             "start"_a, "prot"_a, "length"_a = 0, "data"_a = py::bytes(),
             "Map a range of memory to the virtual address space of the "
             "simulator.\n\nMaps data from `data` into memory and fills the "
             "rest with zeros if `length` is greater than the size of `data`. "
             "If `length` is 0, the size of `data` will be used as length.")
        .def_property_readonly(
            "memimg", [](PyRaspsim &sim) { return MemImg(sim); },
            "Get a memory image object")
        .def("__str__", &PyRaspsim::str,
             "Get the string representation of the current state of the "
             "simulator");
}

void Raspsim::propagate_x86_exception(byte exception, W32 errorcode,
                                      Waddr virtaddr) {
    PyRaspsim::X86Exception = {exception, errorcode, virtaddr, getContext()};
    longjmp(PyRaspsim::simexit, 1);
}

void Raspsim::handle_syscall_32bit(int semantics) {
    // All 32-bit syscalls are unsupported in simulator-only mode
    PyRaspsim::X86Exception = {"Syscall not supported"};
    longjmp(PyRaspsim::simexit, 1);
}

void Raspsim::handle_syscall_64bit() {
    PyRaspsim::X86Exception = {"Syscall not supported"};
    longjmp(PyRaspsim::simexit, 1);
}