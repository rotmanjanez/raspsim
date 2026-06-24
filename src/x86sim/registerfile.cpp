#include "x86sim/registerfile.hpp"

namespace x86sim {

RegisterRef RegisterFile::operator[](Register reg) noexcept {
  return RegisterRef(*this, reg);
}

word_t RegisterFile::operator[](Register reg) const noexcept {
  return commitarf[static_cast<int>(reg)];
}

XmmRegisterRef RegisterFile::operator[](XmmRegister reg) noexcept {
  return XmmRegisterRef(*this, reg);
}

XmmValue RegisterFile::operator[](XmmRegister reg) const noexcept {
  const int lo = static_cast<int>(reg);
  return {.lo = commitarf[lo], .hi = commitarf[lo + 1]};
}

RegisterRef& RegisterRef::operator=(word_t value) noexcept {
  registers_->commitarf[static_cast<int>(reg_)] = value;
  return *this;
}

RegisterRef::operator word_t() const noexcept {
  return registers_->commitarf[static_cast<int>(reg_)];
}

XmmRegisterRef& XmmRegisterRef::operator=(XmmValue value) noexcept {
  const int lo = static_cast<int>(reg_);
  registers_->commitarf[lo] = value.lo;
  registers_->commitarf[lo + 1] = value.hi;
  return *this;
}

XmmRegisterRef::operator XmmValue() const noexcept {
  const int lo = static_cast<int>(reg_);
  return {.lo = registers_->commitarf[lo], .hi = registers_->commitarf[lo + 1]};
}

} // namespace x86sim
