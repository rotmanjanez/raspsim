#ifndef X86SIM_REGISTERFILE_HPP
#define X86SIM_REGISTERFILE_HPP

#include <cstddef>
#include <cstdint>
#include <format>

namespace x86sim {

using word_t = std::uint64_t;

enum class Register {
  rax = 0,
  rcx = 1,
  rdx = 2,
  rbx = 3,
  rsp = 4,
  rbp = 5,
  rsi = 6,
  rdi = 7,
  r8 = 8,
  r9 = 9,
  r10 = 10,
  r11 = 11,
  r12 = 12,
  r13 = 13,
  r14 = 14,
  r15 = 15,
  rip = 56,
  flags = 57
};

enum class XmmRegister {
  xmm0 = 16,
  xmm1 = 18,
  xmm2 = 20,
  xmm3 = 22,
  xmm4 = 24,
  xmm5 = 26,
  xmm6 = 28,
  xmm7 = 30,
  xmm8 = 32,
  xmm9 = 34,
  xmm10 = 36,
  xmm11 = 38,
  xmm12 = 40,
  xmm13 = 42,
  xmm14 = 44,
  xmm15 = 46
};

struct XmmValue {
  word_t lo = 0;
  word_t hi = 0;
};

class RegisterRef;
class XmmRegisterRef;

struct RegisterFile {
  static constexpr std::size_t kRegisterCount = 64;

  [[nodiscard]] RegisterRef operator[](Register reg) noexcept;
  [[nodiscard]] word_t operator[](Register reg) const noexcept;
  [[nodiscard]] XmmRegisterRef operator[](XmmRegister reg) noexcept;
  [[nodiscard]] XmmValue operator[](XmmRegister reg) const noexcept;

  [[nodiscard]] word_t* data() noexcept { return commitarf; }
  [[nodiscard]] const word_t* data() const noexcept { return commitarf; }

  // Public for the existing simulator internals. External API should prefer operator[].
  word_t commitarf[kRegisterCount]{};
};

class RegisterRef {
public:
  RegisterRef& operator=(word_t value) noexcept;
  [[nodiscard]] operator word_t() const noexcept;

private:
  friend struct RegisterFile;
  constexpr RegisterRef(RegisterFile& registers, Register reg) noexcept : registers_(&registers), reg_(reg) {}

  RegisterFile* registers_;
  Register reg_;
};

class XmmRegisterRef {
public:
  XmmRegisterRef& operator=(XmmValue value) noexcept;
  [[nodiscard]] operator XmmValue() const noexcept;

private:
  friend struct RegisterFile;
  constexpr XmmRegisterRef(RegisterFile& registers, XmmRegister reg) noexcept : registers_(&registers), reg_(reg) {}

  RegisterFile* registers_;
  XmmRegister reg_;
};

} // namespace x86sim

template<>
struct std::formatter<x86sim::RegisterRef> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::RegisterRef& reg, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "0x{:016x}", static_cast<x86sim::word_t>(reg));
  }
};

template<>
struct std::formatter<x86sim::XmmRegisterRef> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::XmmRegisterRef& reg, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "{}", static_cast<x86sim::XmmValue>(reg));
  }
};

template<>
struct std::formatter<x86sim::XmmValue> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::XmmValue& value, FormatContext& ctx) const {
    return std::format_to(ctx.out(), "0x{:016x}{:016x}", value.hi, value.lo);
  }
};

template<>
struct std::formatter<x86sim::Register> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::Register reg, FormatContext& ctx) const {
    using enum x86sim::Register;
    switch (reg) {
    case rax:
      return std::format_to(ctx.out(), "rax");
    case rcx:
      return std::format_to(ctx.out(), "rcx");
    case rdx:
      return std::format_to(ctx.out(), "rdx");
    case rbx:
      return std::format_to(ctx.out(), "rbx");
    case rsp:
      return std::format_to(ctx.out(), "rsp");
    case rbp:
      return std::format_to(ctx.out(), "rbp");
    case rsi:
      return std::format_to(ctx.out(), "rsi");
    case rdi:
      return std::format_to(ctx.out(), "rdi");
    case r8:
      return std::format_to(ctx.out(), "r8");
    case r9:
      return std::format_to(ctx.out(), "r9");
    case r10:
      return std::format_to(ctx.out(), "r10");
    case r11:
      return std::format_to(ctx.out(), "r11");
    case r12:
      return std::format_to(ctx.out(), "r12");
    case r13:
      return std::format_to(ctx.out(), "r13");
    case r14:
      return std::format_to(ctx.out(), "r14");
    case r15:
      return std::format_to(ctx.out(), "r15");
    case rip:
      return std::format_to(ctx.out(), "rip");
    case flags:
      return std::format_to(ctx.out(), "flags");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(reg));
  }
};

template<>
struct std::formatter<x86sim::XmmRegister> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::XmmRegister reg, FormatContext& ctx) const {
    const int index = (static_cast<int>(reg) - static_cast<int>(x86sim::XmmRegister::xmm0)) / 2;
    if (index >= 0 && index <= 15 && static_cast<int>(reg) == static_cast<int>(x86sim::XmmRegister::xmm0) + index * 2)
      return std::format_to(ctx.out(), "xmm{}", index);
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(reg));
  }
};

template<>
struct std::formatter<x86sim::RegisterFile> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(const x86sim::RegisterFile& registers, FormatContext& ctx) const {
    static constexpr x86sim::Register regs[] = {
        x86sim::Register::rax, x86sim::Register::rbx, x86sim::Register::rcx, x86sim::Register::rdx,
        x86sim::Register::rsp, x86sim::Register::rbp, x86sim::Register::rsi, x86sim::Register::rdi,
        x86sim::Register::r8,  x86sim::Register::r9,  x86sim::Register::r10, x86sim::Register::r11,
        x86sim::Register::r12, x86sim::Register::r13, x86sim::Register::r14, x86sim::Register::r15,
        x86sim::Register::rip, x86sim::Register::flags,
    };

    auto out = std::format_to(ctx.out(), "{{");
    const char* sep = "";
    for (x86sim::Register reg : regs) {
      out = std::format_to(out, "{}{}=0x{:016x}", sep, reg, registers[reg]);
      sep = " ";
    }
    return std::format_to(out, "}}");
  }
};

#endif
