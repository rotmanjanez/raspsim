// -*- c++ -*-
//
// PTLsim: Cycle Accurate x86-64 Simulator
// Generic interface required by simulator
//
// Copyright 2020-2020 Alexis Engelke <engelke@in.tum.de>
//

#ifndef _PTLSIM_API_H_
#define _PTLSIM_API_H_

#include "x86sim/x86sim.hpp"

#include "globals.h"
#include "ptlhwdef.h"

namespace x86sim {

#define contextcount (1)
#define MAX_CONTEXTS 1

static const Waddr INVALID_PHYSADDR = 0;

extern W64 loadphys(Waddr addr);
extern W64 storemask(Waddr addr, W64 data, byte bytemask);

inline CpuidResult default_cpuid(W32 func, W32, int vcpuid) {
  static const char cpuid_vendor[12 + 1] = "GenuineIntel";
  static const char cpuid_description[48 + 1] = "Intel(R) Xeon(TM) CPU 2.00 GHz                  ";

  auto word = [](const char* str) -> W32 {
    return (static_cast<W32>(static_cast<unsigned char>(str[0])) << 0) |
           (static_cast<W32>(static_cast<unsigned char>(str[1])) << 8) |
           (static_cast<W32>(static_cast<unsigned char>(str[2])) << 16) |
           (static_cast<W32>(static_cast<unsigned char>(str[3])) << 24);
  };

  constexpr W32 x86_feature_fpu = (1 << 0);   // Onboard FPU
  constexpr W32 x86_feature_pse = (1 << 3);   // Page Size Extensions
  constexpr W32 x86_feature_tsc = (1 << 4);   // Time Stamp Counter
  constexpr W32 x86_feature_msr = (1 << 5);   // Model-Specific Registers, RDMSR, WRMSR
  constexpr W32 x86_feature_pae = (1 << 6);   // Physical Address Extensions
  constexpr W32 x86_feature_mce = (1 << 7);   // Machine Check Architecture
  constexpr W32 x86_feature_cx8 = (1 << 8);   // CMPXCHG8 instruction
  constexpr W32 x86_feature_apic = (1 << 9);  // Onboard APIC
  constexpr W32 x86_feature_mca = (1 << 14);  // Machine Check Architecture
  constexpr W32 x86_feature_cmov = (1 << 15); // CMOV instruction
  constexpr W32 x86_feature_pat = (1 << 16);  // Page Attribute Table
  constexpr W32 x86_feature_pse36 = (1 << 17);
  constexpr W32 x86_feature_pn = (1 << 18);    // Processor serial number
  constexpr W32 x86_feature_clfl = (1 << 19);  // Supports the CLFLUSH instruction
  constexpr W32 x86_feature_nx = (1 << 20);    // No-Execute page attribute
  constexpr W32 x86_feature_acpi = (1 << 22);  // ACPI via MSR
  constexpr W32 x86_feature_mmx = (1 << 23);   // Multimedia Extensions
  constexpr W32 x86_feature_fxsr = (1 << 24);  // FXSAVE and FXRSTOR
  constexpr W32 x86_feature_xmm = (1 << 25);   // Streaming SIMD Extensions
  constexpr W32 x86_feature_xmm2 = (1 << 26);  // Streaming SIMD Extensions-2
  constexpr W32 x86_feature_snoop = (1 << 27); // CPU self snoop
  constexpr W32 x86_feature_ht = (1 << 28);    // Hyper-Threading

  constexpr W32 ptlsim_x86_feature =
      x86_feature_fpu | x86_feature_pse | x86_feature_tsc | x86_feature_msr | x86_feature_pae | x86_feature_mce |
      x86_feature_cx8 | x86_feature_apic | x86_feature_mca | x86_feature_cmov | x86_feature_pat | x86_feature_pse36 |
      x86_feature_pn | x86_feature_clfl | x86_feature_nx | x86_feature_acpi | x86_feature_mmx | x86_feature_fxsr |
      x86_feature_xmm | x86_feature_xmm2 | x86_feature_snoop | x86_feature_ht;

  constexpr W32 x86_ext_feature_xmm3 = (1 << 0);  // Streaming SIMD Extensions-3
  constexpr W32 x86_ext_feature_cx16 = (1 << 13); // CMPXCHG16B
  constexpr W32 ptlsim_x86_ext_feature = x86_ext_feature_xmm3 | x86_ext_feature_cx16;

  constexpr W32 x86_vendor_feature_fxsr_opt = (1 << 25); // FXSR optimizations
  constexpr W32 x86_vendor_feature_lm = (1 << 29);       // Long Mode (x86-64)
  constexpr W32 ptlsim_x86_vendor_feature =
      x86_vendor_feature_fxsr_opt | x86_vendor_feature_lm | (ptlsim_x86_feature & 0x1ffffff);

  constexpr W32 x86_vendor_ext_feature_lahf_lm = (1 << 0); // LAHF/SAHF in long mode
  constexpr W32 ptlsim_x86_vendor_ext_feature = x86_vendor_ext_feature_lahf_lm;

  constexpr W32 ptlsim_x86_model_info = (0 << 0) | (0 << 4) | (15 << 8) | (0 << 12) | (0 << 16) | (0 << 20) | (0 << 24);
  constexpr W32 ptlsim_x86_misc_info = (0 << 0) | (8 << 8) | (0 << 16) | (0 << 24);

  switch (func) {
  case 0:
    return {2, word(cpuid_vendor), word(cpuid_vendor + 8), word(cpuid_vendor + 4)};
  case 1:
    return {ptlsim_x86_model_info, ptlsim_x86_misc_info | (static_cast<W32>(vcpuid) << 24), ptlsim_x86_ext_feature,
            ptlsim_x86_feature};
  case 2:
    return {0, 0, 0, 0};
  case 0x80000000:
    return {4, word(cpuid_vendor), word(cpuid_vendor + 8), word(cpuid_vendor + 4)};
  case 0x80000001:
    return {ptlsim_x86_model_info, 0, ptlsim_x86_vendor_ext_feature, ptlsim_x86_vendor_feature};
  case 0x80000002 ... 0x80000004: {
    const char* cpudesc = &cpuid_description[(func - 0x80000002) * 16];
    return {word(cpudesc), word(cpudesc + 4), word(cpudesc + 8), word(cpudesc + 12)};
  }
  default:
    return {0, 0, 0, 0};
  }
}

//
// System calls
//
enum { SYSCALL_SEMANTICS_INT80, SYSCALL_SEMANTICS_SYSCALL, SYSCALL_SEMANTICS_SYSENTER };

void handle_syscall_32bit(Context& context, int semantics);

// x86-64 mode has only one type of system call (the syscall instruction)
void handle_syscall_64bit(Context& context);

CpuidResult handle_cpuid(Context& context, W32 func, W32 subfunc);

// Used to determine whether to exit emulation.
extern bool requested_switch_to_native;


} // namespace x86sim

#endif // _PTLSIM_API_H_
