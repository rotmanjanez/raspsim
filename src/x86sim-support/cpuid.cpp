#include "x86sim-support/defaults.hpp"

namespace x86sim::defaults {

CpuidResult default_cpuid(CpuidRequest request, int vcpuid) noexcept {
  static const char cpuid_vendor[12 + 1] = "GenuineIntel";
  static const char cpuid_description[48 + 1] = "Intel(R) Xeon(TM) CPU 2.00 GHz                  ";

  auto word = [](const char* str) -> std::uint32_t {
    return (static_cast<std::uint32_t>(static_cast<unsigned char>(str[0])) << 0) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(str[1])) << 8) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(str[2])) << 16) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(str[3])) << 24);
  };

  constexpr std::uint32_t x86_feature_fpu = (1 << 0);
  constexpr std::uint32_t x86_feature_pse = (1 << 3);
  constexpr std::uint32_t x86_feature_tsc = (1 << 4);
  constexpr std::uint32_t x86_feature_msr = (1 << 5);
  constexpr std::uint32_t x86_feature_pae = (1 << 6);
  constexpr std::uint32_t x86_feature_mce = (1 << 7);
  constexpr std::uint32_t x86_feature_cx8 = (1 << 8);
  constexpr std::uint32_t x86_feature_apic = (1 << 9);
  constexpr std::uint32_t x86_feature_mca = (1 << 14);
  constexpr std::uint32_t x86_feature_cmov = (1 << 15);
  constexpr std::uint32_t x86_feature_pat = (1 << 16);
  constexpr std::uint32_t x86_feature_pse36 = (1 << 17);
  constexpr std::uint32_t x86_feature_pn = (1 << 18);
  constexpr std::uint32_t x86_feature_clfl = (1 << 19);
  constexpr std::uint32_t x86_feature_nx = (1 << 20);
  constexpr std::uint32_t x86_feature_acpi = (1 << 22);
  constexpr std::uint32_t x86_feature_mmx = (1 << 23);
  constexpr std::uint32_t x86_feature_fxsr = (1 << 24);
  constexpr std::uint32_t x86_feature_xmm = (1 << 25);
  constexpr std::uint32_t x86_feature_xmm2 = (1 << 26);
  constexpr std::uint32_t x86_feature_snoop = (1 << 27);
  constexpr std::uint32_t x86_feature_ht = (1 << 28);

  constexpr std::uint32_t ptlsim_x86_feature =
      x86_feature_fpu | x86_feature_pse | x86_feature_tsc | x86_feature_msr | x86_feature_pae | x86_feature_mce |
      x86_feature_cx8 | x86_feature_apic | x86_feature_mca | x86_feature_cmov | x86_feature_pat |
      x86_feature_pse36 | x86_feature_pn | x86_feature_clfl | x86_feature_nx | x86_feature_acpi | x86_feature_mmx |
      x86_feature_fxsr | x86_feature_xmm | x86_feature_xmm2 | x86_feature_snoop | x86_feature_ht;

  constexpr std::uint32_t x86_ext_feature_xmm3 = (1 << 0);
  constexpr std::uint32_t x86_ext_feature_cx16 = (1 << 13);
  constexpr std::uint32_t ptlsim_x86_ext_feature = x86_ext_feature_xmm3 | x86_ext_feature_cx16;

  constexpr std::uint32_t x86_vendor_feature_fxsr_opt = (1 << 25);
  constexpr std::uint32_t x86_vendor_feature_lm = (1 << 29);
  constexpr std::uint32_t ptlsim_x86_vendor_feature =
      x86_vendor_feature_fxsr_opt | x86_vendor_feature_lm | (ptlsim_x86_feature & 0x1ffffff);

  constexpr std::uint32_t x86_vendor_ext_feature_lahf_lm = (1 << 0);
  constexpr std::uint32_t ptlsim_x86_vendor_ext_feature = x86_vendor_ext_feature_lahf_lm;

  constexpr std::uint32_t ptlsim_x86_model_info =
      (0 << 0) | (0 << 4) | (15 << 8) | (0 << 12) | (0 << 16) | (0 << 20) | (0 << 24);
  constexpr std::uint32_t ptlsim_x86_misc_info = (0 << 0) | (8 << 8) | (0 << 16) | (0 << 24);

  switch (request.function) {
  case 0:
    return {2, word(cpuid_vendor), word(cpuid_vendor + 8), word(cpuid_vendor + 4)};
  case 1:
    return {ptlsim_x86_model_info, ptlsim_x86_misc_info | (static_cast<std::uint32_t>(vcpuid) << 24),
            ptlsim_x86_ext_feature, ptlsim_x86_feature};
  case 2:
    return {0, 0, 0, 0};
  case 0x80000000:
    return {4, word(cpuid_vendor), word(cpuid_vendor + 8), word(cpuid_vendor + 4)};
  case 0x80000001:
    return {ptlsim_x86_model_info, 0, ptlsim_x86_vendor_ext_feature, ptlsim_x86_vendor_feature};
  case 0x80000002:
  case 0x80000003:
  case 0x80000004: {
    const char* cpudesc = &cpuid_description[(request.function - 0x80000002) * 16];
    return {word(cpudesc), word(cpudesc + 4), word(cpudesc + 8), word(cpudesc + 12)};
  }
  default:
    return {0, 0, 0, 0};
  }
}

} // namespace x86sim::defaults
