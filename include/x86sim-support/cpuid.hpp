#ifndef X86SIM_SUPPORT_CPUID_HPP
#define X86SIM_SUPPORT_CPUID_HPP

#include "x86sim/x86sim.hpp"

namespace x86sim::defaults {

[[nodiscard]] CpuidResult default_cpuid(CpuidRequest request, int vcpuid = 0) noexcept;

} // namespace x86sim::defaults

#endif
