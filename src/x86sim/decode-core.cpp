//
// PTLsim: Cycle Accurate x86-64 Simulator
// Decoder for x86 and x86-64 to PTL transops
//
// Copyright 1999-2008 Matt T. Yourst <yourst@yourst.com>
//

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <format>
#include <limits>
#include <string>
#include "globals.h"
#include "ptlhwdef.h"
#include "ptlsim.h"
#include "decode.h"
#include "stats.h"
#include "typedefs.h"
#include "x86sim/addrspace.hpp"
#include "x86sim/logging.hpp"
#include "x86sim/x86sim.hpp"


namespace x86sim {

FILE* bbcache_dump_file = nullptr;

//
// Calling convention:
// rip = return RIP after insn
// sr0 = RIP of insn
// sr1 = RIP after insn
// sr2 = argument
//

const assist_func_t assistid_to_func[ASSIST_COUNT] = {
    // Forced assists based on decode context
    assist_invalid_opcode,
    assist_exec_page_fault,
    assist_gp_fault,
    // x87
    assist_x87_fist,
    assist_x87_fldcw,
    assist_x87_fprem,
    assist_x87_fyl2xp1,
    assist_x87_fsqrt,
    assist_x87_fsincos,
    assist_x87_frndint,
    assist_x87_fscale,
    assist_x87_fsin,
    assist_x87_fcos,
    assist_x87_fxam,
    assist_x87_f2xm1,
    assist_x87_fyl2x,
    assist_x87_fptan,
    assist_x87_fpatan,
    assist_x87_fxtract,
    assist_x87_fprem1,
    assist_x87_fld80,
    assist_x87_fstp80,
    assist_x87_fsave,
    assist_x87_frstor,
    assist_x87_finit,
    assist_x87_fclex,
    // SSE save/restore
    assist_ldmxcsr,
    assist_fxsave,
    assist_fxrstor,
    // Interrupts, system calls, etc.
    assist_int,
    assist_syscall,
    assist_hypercall,
    assist_ptlcall,
    assist_sysenter,
    assist_iret16,
    assist_iret32,
    assist_iret64,
    // Control register updates
    assist_cpuid,
    assist_rdtsc,
    assist_cld,
    assist_std,
    assist_popf,
    assist_write_segreg,
    assist_wrmsr,
    assist_rdmsr,
    assist_write_cr0,
    assist_write_cr2,
    assist_write_cr3,
    assist_write_cr4,
    assist_write_debug_reg,
    // I/O and legacy
    assist_ioport_in,
    assist_ioport_out,
};

int assist_index(assist_func_t assist) {
  foreach (i, ASSIST_COUNT) {
    if (assistid_to_func[i] == assist) {
      return i;
    }
  }

  return -1;
}

const char* assist_name(assist_func_t assist) {
  foreach (i, ASSIST_COUNT) {
    if (assistid_to_func[i] == assist) {
      return assist_names[i];
    }
  }

  return "unknown";
}

void update_assist_stats(assist_func_t assist) {
  int idx = assist_index(assist);
  assert(inrange(idx, 0, ASSIST_COUNT - 1));
  stats.external.assists[idx]++;
}

void split_unaligned(const TransOp& transop, TransOpBuffer& buf) {
  assert(transop.unaligned);

  bool ld = isload(transop.opcode);
  bool st = isstore(transop.opcode);

  assert(ld | st);

  buf.reset();

  int idx;

  idx = buf.put();
  TransOp& ag = buf.uops[idx];
  ag = transop;
  ag.opcode = OP_add;
  ag.size = 3;
  ag.cond = 0;
  ag.eom = 0;
  ag.internal = 0;
  ag.unaligned = 0;
  ag.rd = REG_temp9;
  ag.rc = REG_zero;
  buf.synthops[idx] = get_synthcode_for_uop(OP_add, 3, 0, 0, 0, 0, 0);

  idx = buf.put();
  TransOp& lo = buf.uops[idx];
  lo = transop;
  lo.ra = REG_temp9;
  lo.rb = REG_imm;
  lo.rbimm = 0;
  lo.cond = LDST_ALIGN_LO;
  lo.unaligned = 0;
  lo.som = 0;
  lo.eom = 0;
  buf.synthops[idx] = {}; // loads and stores are not synthesized

  idx = buf.put();
  TransOp& hi = buf.uops[idx];
  hi = transop;
  hi.ra = REG_temp9;
  hi.rb = REG_imm;
  hi.rbimm = 0;
  hi.cond = LDST_ALIGN_HI;
  hi.unaligned = 0;
  hi.som = 0;
  buf.synthops[idx] = {}; // loads and stores are not synthesized

  if (ld) {
    // ld rd = [ra+rb]        =>   ld.lo rd = [rt]           and    ld.hi rd = [rt],rd
    lo.rd = REG_temp4;
    lo.size = 3; // always load 64-bit word
    hi.rb = REG_temp4;
  } else {
    //
    // For stores, expand     st sfrd = [ra+rb],rc    =>   st.lo sfrd1 = [rt],rc    and    st.hi sfrd2 = [rt],rc
    //
    // Make sure high part issues first in program order, so if there is
    // a page fault on the high page it overlaps, this will be triggered
    // before the low part overwrites memory.
    //
    lo.cond = LDST_ALIGN_HI;
    hi.cond = LDST_ALIGN_LO;
  }
}

// clang-format off
static const W16 prefix_map_x86_64[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, PFX_ES, 0, 0, 0, 0, 0, 0, 0, PFX_CS, 0,
  0, 0, 0, 0, 0, 0, PFX_SS, 0, 0, 0, 0, 0, 0, 0, PFX_DS, 0,
  PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX, PFX_REX,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, PFX_FS, PFX_GS, PFX_DATA, PFX_ADDR, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, PFX_FWAIT, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  PFX_LOCK, 0, PFX_REPNZ, PFX_REPZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

static const W16 prefix_map_x86[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, PFX_ES, 0, 0, 0, 0, 0, 0, 0, PFX_CS, 0,
  0, 0, 0, 0, 0, 0, PFX_SS, 0, 0, 0, 0, 0, 0, 0, PFX_DS, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, PFX_FS, PFX_GS, PFX_DATA, PFX_ADDR, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, PFX_FWAIT, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  PFX_LOCK, 0, PFX_REPNZ, PFX_REPZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

const char* prefix_names[PFX_count] = {"repz", "repnz", "lock", "cs", "ss", "ds", "es", "fs", "gs", "datasz", "addrsz", "rex", "fwait"};

const char* uniform_arch_reg_names[APR_COUNT] = {
  // 64-bit
  "rax", "rcx", "rdx", "rbx", "rsp", "rbp", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15",
  // 32-bit
  "eax", "ecx", "edx", "ebx", "esp", "ebp", "esi", "edi", "r8d", "r9d", "r10d", "r11d", "r12d", "r13d", "r14d", "r15d",
  // 16-bit
  "ax", "cx", "dx", "bx", "sp", "bp", "si", "di", "r8w", "r9w", "r10w", "r11w", "r12w", "r13w", "r14w", "r15w",
  // 8-bit
  "al", "cl", "dl", "bl", "ah", "ch", "dh", "bh",
  // 8-bit with REX:
  "spl", "bpl", "sil", "dil",
  "r8b", "r9b", "r10b", "r11b", "r12b", "r13b", "r14b", "r15b",
  // SSE registers
  "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15",
  // segments
  "es", "cs", "ss", "ds", "fs", "gs",
  // special:
  "rip", "zero",
};

const ArchPseudoRegInfo reginfo[APR_COUNT] = {
  // 64-bit
  {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0},
  // 32-bit
  {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0}, {2, 0},
  // 16-bit
  {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0},
  // 8-bit
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 1}, {0, 1}, {0, 1}, {0, 1},
  // 8-bit with REX, not double-counting the regular 8-bit regs:
  {0, 0}, {0, 0}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
  // SSE registers
  {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0}, {3, 0},
  // segments:
  {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0},
  // special:
  {3, 0}, {3, 0},
};

const byte reg64_to_uniform_reg[16] = { APR_rax, APR_rcx, APR_rdx, APR_rbx, APR_rsp, APR_rbp, APR_rsi, APR_rdi, APR_r8, APR_r9, APR_r10, APR_r11, APR_r12, APR_r13, APR_r14, APR_r15 };
const byte xmmreg_to_uniform_reg[16] = { APR_xmm0, APR_xmm1, APR_xmm2, APR_xmm3, APR_xmm4, APR_xmm5, APR_xmm6, APR_xmm7, APR_xmm8, APR_xmm9, APR_xmm10, APR_xmm11, APR_xmm12, APR_xmm13, APR_xmm14, APR_xmm15 };
const byte reg32_to_uniform_reg[16] = { APR_eax, APR_ecx, APR_edx, APR_ebx, APR_esp, APR_ebp, APR_esi, APR_edi, APR_r8d, APR_r9d, APR_r10d, APR_r11d, APR_r12d, APR_r13d, APR_r14d, APR_r15d };
const byte reg16_to_uniform_reg[16] = { APR_ax, APR_cx, APR_dx, APR_bx, APR_sp, APR_bp, APR_si, APR_di, APR_r8w, APR_r9w, APR_r10w, APR_r11w, APR_r12w, APR_r13w, APR_r14w, APR_r15w };
const byte reg8_to_uniform_reg[8] = { APR_al, APR_cl, APR_dl, APR_bl, APR_ah, APR_ch, APR_dh, APR_bh };
const byte reg8x_to_uniform_reg[16] = { APR_al, APR_cl, APR_dl, APR_bl, APR_spl, APR_bpl, APR_sil, APR_dil, APR_r8b, APR_r9b, APR_r10b, APR_r11b, APR_r12b, APR_r13b, APR_r14b, APR_r15b };
const byte segreg_to_uniform_reg[16] = { APR_es, APR_cs, APR_ss, APR_ds, APR_fs, APR_zero, APR_zero };

const byte arch_pseudo_reg_to_arch_reg[APR_COUNT] = {
  // 64-bit
  REG_rax, REG_rcx, REG_rdx, REG_rbx, REG_rsp, REG_rbp, REG_rsi, REG_rdi, REG_r8, REG_r9, REG_r10, REG_r11, REG_r12, REG_r13, REG_r14, REG_r15,
  // 32-bit
  REG_rax, REG_rcx, REG_rdx, REG_rbx, REG_rsp, REG_rbp, REG_rsi, REG_rdi, REG_r8, REG_r9, REG_r10, REG_r11, REG_r12, REG_r13, REG_r14, REG_r15,
  // 16-bit
  REG_rax, REG_rcx, REG_rdx, REG_rbx, REG_rsp, REG_rbp, REG_rsi, REG_rdi, REG_r8, REG_r9, REG_r10, REG_r11, REG_r12, REG_r13, REG_r14, REG_r15,
  // 8-bit
  REG_rax, REG_rcx, REG_rdx, REG_rbx, REG_rax, REG_rcx, REG_rdx, REG_rbx,
  // 8-bit with REX, not double-counting the regular 8-bit regs:
  REG_rsp, REG_rbp, REG_rsi, REG_rdi,
  REG_r8, REG_r9, REG_r10, REG_r11, REG_r12, REG_r13, REG_r14, REG_r15,
  // SSE registers
  REG_xmml0, REG_xmml1, REG_xmml2, REG_xmml3, REG_xmml4, REG_xmml5, REG_xmml6, REG_xmml7, REG_xmml8, REG_xmml9, REG_xmml10, REG_xmml11, REG_xmml12, REG_xmml13, REG_xmml14, REG_xmml15,
  // segments:
  REG_zero, REG_zero, REG_zero, REG_zero, REG_zero, REG_zero,
  // special:
  REG_rip, REG_zero
};

// For specifying easy to read arrays
#define _ (0)

static const byte onebyte_has_modrm[256] = {
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
  /*       -------------------------------        */
  /* 00 */ 1,1,1,1,_,_,_,_,1,1,1,1,_,_,_,_, /* 00 */
  /* 10 */ 1,1,1,1,_,_,_,_,1,1,1,1,_,_,_,_, /* 10 */
  /* 20 */ 1,1,1,1,_,_,_,_,1,1,1,1,_,_,_,_, /* 20 */
  /* 30 */ 1,1,1,1,_,_,_,_,1,1,1,1,_,_,_,_, /* 30 */
  /* 40 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 40 */
  /* 50 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 50 */
  /* 60 */ _,_,1,1,_,_,_,_,_,1,_,1,_,_,_,_, /* 60 */
  /* 70 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 70 */
  /* 80 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 80 */
  /* 90 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 90 */
  /* a0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* a0 */
  /* b0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* b0 */
  /* c0 */ 1,1,_,_,1,1,1,1,_,_,_,_,_,_,_,_, /* c0 */
  /* d0 */ 1,1,1,1,_,_,_,_,1,1,1,1,1,1,1,1, /* d0 */
  /* e0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* e0 */
  /* f0 */ _,_,_,_,_,_,1,1,_,_,_,_,_,_,1,1  /* f0 */
  /*       -------------------------------        */
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
};

static const byte twobyte_has_modrm[256] = {
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
  /*       -------------------------------        */
  /* 00 */ 1,1,1,1,_,_,_,_,_,_,_,_,_,1,_,1, /* 0f */
  /* 10 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 1f */
  /* 20 */ 1,1,1,1,1,_,1,_,1,1,1,1,1,1,1,1, /* 2f */
  /* 30 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 3f */
  /* 40 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 4f */
  /* 50 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 5f */
  /* 60 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 6f */
  /* 70 */ 1,1,1,1,1,1,1,_,_,_,_,_,1,1,1,1, /* 7f */
  /* 80 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 8f */
  /* 90 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 9f */
  /* a0 */ _,_,_,1,1,1,1,1,_,_,_,1,1,1,1,1, /* af */
  /* b0 */ 1,1,1,1,1,1,1,1,_,_,1,1,1,1,1,1, /* bf */
  /* c0 */ 1,1,1,1,1,1,1,1,_,_,_,_,_,_,_,_, /* cf */
  /* d0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* df */
  /* e0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* ef */
  /* f0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,_  /* ff */
  /*       -------------------------------        */
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
};

static const byte twobyte_uses_SSE_prefix[256] = {
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
  /*       -------------------------------        */
  /* 00 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 0f */
  /* 10 */ 1,1,1,1,1,1,1,1,_,_,_,_,_,_,_,_, /* 1f */
  /* 20 */ _,_,_,_,_,_,_,_,1,1,1,1,1,1,1,1, /* 2f */
  /* 30 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 3f */
  /* 40 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 4f */
  /* 50 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 5f */
  /* 60 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 6f */
  /* 70 */ 1,1,1,1,1,1,1,_,_,_,_,_,1,1,1,1, /* 7f */
  /* 80 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 8f */
  /* 90 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* 9f */
  /* a0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* af */
  /* b0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /* bf */
  /* c0 */ _,_,1,_,1,1,1,_,_,_,_,_,_,_,_,_, /* cf */
  /* d0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* df */
  /* e0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* ef */
  /* f0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,_  /* ff */
  /*       -------------------------------        */
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
};

//
// This determines if the insn is handled by the
// fast decoder or the complex microcode decoder.
// The expanded x86 opcodes are from 0x000 to 0x1ff,
// i.e. normal ones and those with the 0x0f prefix:
//
static const byte insn_is_simple[512] = {
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
  /*       -------------------------------        */
  /* 00 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,_, /* 0f */
  /* 10 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 1f */
  /* 20 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 2f */
  /* 30 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 3f */
  /* 40 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 4f */
  /* 50 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 5f */
  /* 60 */ _,_,_,1,_,_,_,_,1,1,1,1,_,_,_,_, /* 6f */
  /* 70 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 7f */
  /* 80 */ 1,1,1,1,1,1,_,_,1,1,1,1,_,1,_,1, /* 8f */
  /* 90 */ 1,_,_,_,_,_,_,_,1,1,_,_,_,_,1,1, /* 9f */
  /* a0 */ 1,1,1,1,_,_,_,_,1,1,_,_,_,_,_,_, /* af */
  /* b0 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* bf */
  /* c0 */ 1,1,1,1,_,_,1,1,1,1,_,_,_,_,_,_, /* cf */
  /* d0 */ 1,1,1,1,_,_,_,_,_,_,_,_,_,_,_,_, /* df */
  /* e0 */ _,_,_,_,_,_,_,_,1,1,_,1,_,_,_,_, /* ef */
  /* f0 */ _,_,_,_,_,1,1,1,1,1,_,_,0,0,1,1, /* ff */
  /*100 */ _,_,_,_,_,_,_,_,_,_,_,_,_,1,_,_, /*10f */
  /*110 */ _,_,_,_,_,_,_,_,1,_,_,_,_,_,_,_, /*11f */
  /*120 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*12f */
  /*130 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*13f */
  /*140 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /*14f */
  /*150 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*15f */
  /*160 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*16f */
  /*170 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*17f */
  /*180 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /*18f */
  /*190 */ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /*19f */
  /*1a0 */ _,_,_,1,_,_,_,_,_,_,_,1,_,_,_,1, /*1af */
  /*1b0 */ _,_,_,1,_,_,1,1,_,_,1,1,1,1,1,1, /*1bf */
  /*1c0 */ _,_,_,_,_,_,_,_,1,1,1,1,1,1,1,1, /*1cf */
  /*1d0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*1df */
  /*1e0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*1ef */
  /*1f0 */ _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, /*1ff */
  /*       -------------------------------        */
  /*       0 1 2 3 4 5 6 7 8 9 a b c d e f        */
};
// clang-format on
#undef _

static int transop_histogram[MAX_TRANSOPS_PER_USER_INSN + 1];

void TraceDecoder::reset() {
  byteoffset = 0;
  transbufcount = 0;
  pfec = 0;
  faultaddr = 0;
  prefixes = 0;
  rex = 0;
  modrm = 0;
  user_insn_count = 0;
  // Note(AE): commit d278b597 changed this to zero, previously this was always
  // reset to 1 from the very beginning. The commit mentions that this was done
  // to ensure that every basic block updates flags at least once, but does not
  // give any reason/indication of *why* this was done. I'm unsure whether
  // anything depends on this behavior now, but at least theoretically, there
  // should be no harm in setting this to one.
  last_flags_update_was_atomic = 1;
  invalid = 0;
  some_insns_complex = 0;
  used_microcode_assist = 0;
  end_of_block = 0;
  first_uop_in_insn = 1;
  split_basic_block_at_locks_and_fences = 0;
  split_invalid_basic_blocks = 0;
  no_partial_flag_updates_per_insn = 0;
  fast_length_decode_only = 0;
  join_with_prev_insn = 0;
  outcome = DECODE_OUTCOME_OK;
  stop_at_rip = std::numeric_limits<W64>::max();
  stop_at_user_insns = std::numeric_limits<W64>::max();
  machine = nullptr;
}

TraceDecoder::TraceDecoder(const RIPVirtPhys& rvp) {
  reset();
  bb.reset(rvp);
  rip = rvp;
  ripstart = rvp;
  use64 = rvp.use64;
  kernel = rvp.kernel;
  dirflag = rvp.df;
}

TraceDecoder::TraceDecoder(Waddr rip, bool use64, bool kernel, bool df) {
  reset();
  bb.reset();
  setzero(bb.rip);
  bb.rip.rip = rip;
  bb.rip.use64 = use64;
  bb.rip.kernel = kernel;
  bb.rip.df = df;
  this->rip = rip;
  this->ripstart = rip;
  this->use64 = use64;
  this->kernel = kernel;
  this->dirflag = df;
}

TraceDecoder::TraceDecoder(Context& ctx, Waddr rip) {
  reset();
  use64 = ctx.use64;
  kernel = 0;
  dirflag = ((ctx.internal_eflags & FLAG_DF) != 0);

  bb.reset();
  setzero(bb.rip);
  bb.rip.rip = rip;
  bb.rip.use64 = use64;
  bb.rip.kernel = kernel;
  bb.rip.df = dirflag;

  this->rip = rip;
  ripstart = rip;
}

void TraceDecoder::put(const TransOp& transop) {
  assert(transbufcount < MAX_TRANSOPS_PER_USER_INSN);
  transbuf[transbufcount++] = transop;

  if (transop.setflags)
    last_flags_update_was_atomic = (transop.setflags == 0x7);
}

bool TraceDecoder::flush() {
  if unlikely (!transbufcount) {
    return true;
  }

  //
  // Can we fit all the uops in this x86 insn?
  //
  // If not, remove the pending uops and split
  // the BB with an uncond branch.
  //
  // We reserve 2 uops for the splitting branch
  // and a possible collcc if needed.
  //
  bool overflow = (transbufcount >= ((MAX_BB_UOPS - 2) - bb.count));

  if unlikely (overflow) {
    logging::println(
        logging::DEBUG,
        "Basic block overflowed (too many uops) during decode of {} (ripstart {}): req {} uops but only have {} free",
        W64(bb.rip), (void*)ripstart, transbufcount, ((MAX_BB_UOPS - 2) - bb.count));
    assert(!first_insn_in_bb());
    transbufcount = 0;
    // This will set join_with_prev_insn and fill in the (collcc, bru) ops
    split_before();
  }

  if unlikely (join_with_prev_insn) {
    //
    // Reopen the previous instruction
    //
    assert(bb.count > 0);
    int i = bb.count - 1;
    while (i >= 0) {
      TransOp& prevop = bb.transops[i];
      prevop.eom = 0;
      prevop.final_flags_in_insn = 0;
      i--;
      if unlikely (prevop.som) {
        rip = ripstart;
        ripstart -= prevop.bytes;
        break;
      }
    }
  }

  int bytes = (rip - ripstart);
  assert(bytes <= 15);

  TransOp& first = transbuf[0];
  TransOp& last = transbuf[transbufcount - 1];
  first.som = (!join_with_prev_insn);
  last.eom = 1;

  bool contains_branch = isbranch(last.opcode);

  W8s final_archreg_writer[ARCHREG_COUNT];
  memset(final_archreg_writer, 0xff, sizeof(final_archreg_writer));
  W8s final_flags_writer = -1;

  byte flag_sets_set = 0;
  foreach (i, transbufcount) {
    TransOp& transop = transbuf[i];
    if likely (transop.rd < ARCHREG_COUNT) {
      final_archreg_writer[transop.rd] = i;
    }
    bool sets_all_flags = ((transop.setflags == 7) && (!transop.nouserflags));
    if unlikely (sets_all_flags)
      final_flags_writer = i;
    if likely (!transop.nouserflags)
      flag_sets_set |= transop.setflags;
  }

  if unlikely (no_partial_flag_updates_per_insn) {
    if unlikely ((flag_sets_set != 0) && (flag_sets_set != (SETFLAG_ZF | SETFLAG_CF | SETFLAG_OF))) {
      logging::println(logging::ERROR, "Invalid partial flag sets at rip {} (flag sets {})", (void*)ripstart,
                       flag_sets_set);
      assert(false);
    }
  }

  foreach (i, transbufcount) {
    TransOp& transop = transbuf[i];
    if unlikely (bb.count >= MAX_BB_UOPS) {
      logging::println(logging::ERROR,
                       "ERROR: Too many transops ({}) in basic block {} (current RIP: {}) (max {} allowed)", bb.count,
                       W64(bb.rip), (void*)ripstart, MAX_BB_UOPS);
      logging::println(logging::ERROR, "{}", bb);
      assert(bb.count < MAX_BB_UOPS);
    }

    bool ld = isload(transop.opcode);
    bool st = isstore(transop.opcode);
    bool br = isbranch(transop.opcode);
    if unlikely (br) {
      switch (transop.opcode) {
      case OP_br:
        bb.type = BB_TYPE_COND;
        break;
      case OP_bru:
        bb.type = BB_TYPE_UNCOND;
        break;
      case OP_jmp:
        bb.type = BB_TYPE_INDIR;
        break;
      case OP_brp:
        bb.type = BB_TYPE_ASSIST;
        break;
      default:
        assert(false);
        break;
      }
      bb.call = ((transop.extshift & BRANCH_HINT_PUSH_RAS) != 0);
      bb.ret = ((transop.extshift & BRANCH_HINT_POP_RAS) != 0);
      bb.rip_taken = transop.riptaken;
      bb.rip_not_taken = transop.ripseq;
    }

    transop.unaligned = 0;
    transop.bytes = bytes;
    transop.bbindex = bb.count;
    transop.is_x87 = is_x87;
    transop.is_sse = is_sse;
    transop.final_insn_in_bb = contains_branch;
    transop.final_arch_in_insn = (transop.rd < ARCHREG_COUNT) && (final_archreg_writer[transop.rd] == i);
    transop.final_flags_in_insn = (final_flags_writer == i);
    transop.any_flags_in_insn = (flag_sets_set != 0);

    bb.x87 |= is_x87;
    bb.sse |= is_sse;

    bb.transops[bb.count++] = transop;
    if (ld | st)
      bb.memcount++;
    if (st)
      bb.storecount++;
    bb.tagcount++;

    if (transop.rd < ARCHREG_COUNT)
      setbit(bb.usedregs, transop.rd);
    if (transop.ra < ARCHREG_COUNT)
      setbit(bb.usedregs, transop.ra);
    if (transop.rb < ARCHREG_COUNT)
      setbit(bb.usedregs, transop.rb);
    if (transop.rc < ARCHREG_COUNT)
      setbit(bb.usedregs, transop.rc);
  }

  stats.decoder.throughput.uops += transbufcount;

  if (!join_with_prev_insn) {
    bb.user_insn_count++;
    bb.bytes += bytes;
    stats.decoder.throughput.x86_insns++;
    stats.decoder.throughput.bytes += bytes;
  }

  transbufcount = 0;

  return (!overflow);
}

} // namespace x86sim

auto std::formatter<x86sim::DecodedOperand>::format(const x86sim::DecodedOperand& op, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  switch (op.type) {
  case OPTYPE_REG:
    out = std::format_to(out, "{}", uniform_arch_reg_names[op.reg.reg]);
    break;
  case OPTYPE_IMM:
    out = std::format_to(out, "{:016x}", op.imm.imm);
    break;
  case OPTYPE_MEM:
    out = std::format_to(out, "mem{} [{} + {}*{} + {:016x}{}]", (1 << op.mem.size),
                         uniform_arch_reg_names[op.mem.basereg], uniform_arch_reg_names[op.mem.indexreg],
                         (1 << op.mem.scale), op.mem.offset, (op.mem.riprel ? " riprel" : ""));
    break;
  default:
    break;
  }
  return out;
}

namespace x86sim {

bool DecodedOperand::gform_ext(TraceDecoder& state, int bytemode, int regfield, bool def64, bool in_rex_base) {
  int add = ((in_rex_base) ? state.rex.extbase : state.rex.extreg) ? 8 : 0;

  this->type = OPTYPE_REG;
  switch (bytemode) {
  case b_mode:
    this->reg.reg = (state.rex) ? reg8x_to_uniform_reg[regfield + add] : reg8_to_uniform_reg[regfield];
    break;
  case w_mode:
    this->reg.reg = reg16_to_uniform_reg[regfield + add];
    break;
  case d_mode:
    this->reg.reg = reg32_to_uniform_reg[regfield + add];
    break;
  case q_mode:
    this->reg.reg = reg64_to_uniform_reg[regfield + add];
    break;
  case v_mode:
  case dq_mode:
    this->reg.reg = (state.rex.mode64 | (def64 & (!state.opsize_prefix))) ? reg64_to_uniform_reg[regfield + add]
                    : ((!state.opsize_prefix) | (bytemode == dq_mode))    ? reg32_to_uniform_reg[regfield + add]
                                                                          : reg16_to_uniform_reg[regfield + add];
    break;
  case x_mode:
    this->reg.reg = xmmreg_to_uniform_reg[regfield + add];
    break;
  default:
    return false;
  }

  return true;
}

bool DecodedOperand::gform(TraceDecoder& state, int bytemode) {
  return gform_ext(state, bytemode, state.modrm.reg);
}

bool DecodedOperand::iform(TraceDecoder& state, int bytemode) {
  this->type = OPTYPE_IMM;
  this->imm.imm = 0;

  switch (bytemode) {
  case b_mode:
    this->imm.imm = (W8s)state.fetch1();
    break;
  case q_mode:
    this->imm.imm = (W64s)state.fetch8();
    break;
  case v_mode:
    // NOTE: Even if rex.mode64 is specified, immediates are never longer than 32 bits (except for mov):
    if (state.rex.mode64 | (!state.opsize_prefix)) {
      this->imm.imm = (W32s)state.fetch4();
    } else {
      this->imm.imm = (W16s)state.fetch2();
    }
    break;
  case w_mode:
    this->imm.imm = (W16s)state.fetch2();
    break;
  default:
    return false;
  }

  return true;
}

bool DecodedOperand::iform64(TraceDecoder& state, int bytemode) {
  this->type = OPTYPE_IMM;
  this->imm.imm = 0;

  switch (bytemode) {
  case b_mode:
    this->imm.imm = (W8s)state.fetch1();
    break;
  case q_mode:
    this->imm.imm = (W64s)state.fetch8();
    break;
  case v_mode:
    if (state.rex.mode64) {
      this->imm.imm = (W64s)state.fetch8();
    } else if (state.opsize_prefix) {
      this->imm.imm = (W16s)state.fetch2();
    } else {
      this->imm.imm = (W32s)state.fetch4();
    }
    break;
  case w_mode:
    this->imm.imm = (W16s)state.fetch2();
    break;
  case d_mode:
    this->imm.imm = (W32s)state.fetch4();
    break;
  default:
    return false;
  }
  return true;
}

bool DecodedOperand::eform(TraceDecoder& state, int bytemode) {
  if (state.modrm.mod == 3) {
    return gform_ext(state, bytemode, state.modrm.rm, false, true);
  }

  type = OPTYPE_MEM;
  mem.basereg = APR_zero;
  mem.indexreg = APR_zero;
  mem.offset = 0;
  mem.scale = 0;
  mem.riprel = 0;
  mem.size = 0;

  const int mod_and_rexextbase_and_rm_to_basereg_x86_64[4][2][8] = {
      {
          // mod = 00
          {APR_rax, APR_rcx, APR_rdx, APR_rbx, -1, APR_rip, APR_rsi, APR_rdi}, // rex.extbase = 0
          {APR_r8, APR_r9, APR_r10, APR_r11, -1, APR_rip, APR_r14, APR_r15},   // rex.extbase = 1
      },
      {
          // mod = 01
          {APR_rax, APR_rcx, APR_rdx, APR_rbx, -1, APR_rbp, APR_rsi, APR_rdi}, // rex.extbase = 0
          {APR_r8, APR_r9, APR_r10, APR_r11, -1, APR_r13, APR_r14, APR_r15},   // rex.extbase = 1
      },
      {
          // mod = 10
          {APR_rax, APR_rcx, APR_rdx, APR_rbx, -1, APR_rbp, APR_rsi, APR_rdi}, // rex.extbase = 0
          {APR_r8, APR_r9, APR_r10, APR_r11, -1, APR_r13, APR_r14, APR_r15},   // rex.extbase = 1
      },
      {
          // mod = 11: not possible since this is g-form
          {-1, -1, -1, -1, -1, -1, -1, -1},
          {-1, -1, -1, -1, -1, -1, -1, -1},
      }};

  const int mod_and_rm_to_basereg_x86[4][8] = {
      {APR_eax, APR_ecx, APR_edx, APR_ebx, -1, APR_zero, APR_esi, APR_edi},
      {APR_eax, APR_ecx, APR_edx, APR_ebx, -1, APR_ebp, APR_esi, APR_edi},
      {APR_eax, APR_ecx, APR_edx, APR_ebx, -1, APR_ebp, APR_esi, APR_edi},
      {-1, -1, -1, -1, -1, -1, -1, -1}, // mod = 11: not possible since this is g-form
  };

  mem.basereg = (state.use64)
                    ? mod_and_rexextbase_and_rm_to_basereg_x86_64[state.modrm.mod][state.rex.extbase][state.modrm.rm]
                    : mod_and_rm_to_basereg_x86[state.modrm.mod][state.modrm.rm];

  SIBByte sib;
  if (state.modrm.rm == 4) {
    sib = SIBByte(state.fetch1());
  }

  const byte mod_and_rm_to_immsize[4][8] = {
      {0, 0, 0, 0, 0, 4, 0, 0},
      {1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4},
      {0, 0, 0, 0, 0, 0, 0, 0},
  };

  byte immsize = mod_and_rm_to_immsize[state.modrm.mod][state.modrm.rm];
  mem.offset = (immsize) ? signext32((W32s)state.fetch(immsize), immsize * 8) : 0;
  mem.riprel = (mem.basereg == APR_rip);

  if (mem.basereg < 0) {
    // Have sib
    const int rexextbase_and_base_to_basereg[2][8] = {
        {APR_rax, APR_rcx, APR_rdx, APR_rbx, APR_rsp, -1, APR_rsi, APR_rdi}, // rex.extbase = 0
        {APR_r8, APR_r9, APR_r10, APR_r11, APR_r12, -1, APR_r14, APR_r15},   // rex.extbase = 1
    };

    mem.basereg = rexextbase_and_base_to_basereg[state.rex.extbase][sib.base];
    if (mem.basereg < 0) {
      const int rexextbase_and_mod_to_basereg[2][4] = {
          {APR_zero, APR_rbp, APR_rbp, -1}, // rex.extbase = 0
          {APR_zero, APR_r13, APR_r13, -1}, // rex.extbase = 1
      };

      mem.basereg = rexextbase_and_mod_to_basereg[state.rex.extbase][state.modrm.mod];

      if (!immsize) {
        switch (state.modrm.mod) {
        case 0:
        case 2:
          assert(!immsize);
          mem.offset = (W32s)state.fetch4();
          break;
        case 1:
          assert(!immsize);
          mem.offset = (W8s)state.fetch1();
          break;
        }
      }
    }

    const int rexextindex_and_index_to_indexreg[2][8] = {
        {APR_rax, APR_rcx, APR_rdx, APR_rbx, APR_zero, APR_rbp, APR_rsi, APR_rdi}, // rex.extindex = 0
        {APR_r8, APR_r9, APR_r10, APR_r11, APR_r12, APR_r13, APR_r14, APR_r15},    // rex.extindex = 1
    };

    mem.indexreg = rexextindex_and_index_to_indexreg[state.rex.extindex][sib.index];
    mem.scale = sib.scale;
  }

  switch (bytemode) {
  case b_mode:
    mem.size = 0;
    break;
  case w_mode:
    mem.size = 1;
    break;
  case d_mode:
    mem.size = 2;
    break;
  case q_mode:
    mem.size = 3;
    break;
    // case m_mode: mem.size = (state.use64) ? 3 : 2; break;
  case v_mode:
  case dq_mode:
    mem.size = (state.rex.mode64) ? 3 : ((!state.opsize_prefix) | (bytemode == dq_mode)) ? 2 : 1;
    break; // See table 1.2 (p35) of AMD64 ISA manual
  case x_mode:
    mem.size = 3;
    break;
  default:
    return false;
  }

  return true;
}

bool DecodedOperand::varreg(TraceDecoder& state, int regcode, bool def64) {
  this->type = OPTYPE_REG;

  // push and pop default to 64 bits in 64-bit mode, while all others default to 32 bit mode and need the REX prefix to make them 64-bit:
  // assert(mode_64bit);

  if (def64) {
    // Always a 64-bit operation
    this->reg.reg = reg64_to_uniform_reg[regcode + (state.rex.extbase * 8)];
  } else {
    this->reg.reg = (state.rex.mode64)      ? reg64_to_uniform_reg[regcode + (state.rex.extbase * 8)]
                    : (state.opsize_prefix) ? reg16_to_uniform_reg[regcode + (state.rex.extbase * 8)]
                                            : reg32_to_uniform_reg[regcode + (state.rex.extbase * 8)];
  }

  return true;
}

bool DecodedOperand::varreg_def64(TraceDecoder& state, int regcode) {
  return DecodedOperand::varreg(state, regcode, true);
}

bool DecodedOperand::varreg_def32(TraceDecoder& state, int regcode) {
  return DecodedOperand::varreg(state, regcode, false);
}

void TraceDecoder::immediate(int rdreg, int sizeshift, W64s imm, bool issigned) {
  int totalbits = (sizeshift == 3) ? 64 : (8 * (1 << sizeshift));
  if (totalbits < 64)
    imm = (issigned) ? signext64(imm, totalbits) : bits(imm, 0, totalbits);
  // Only byte and word sized immediates need to be merged with the previous value:
  put(TransOp(OP_mov, rdreg, REG_zero, REG_imm, REG_zero, 3, imm));
}

void TraceDecoder::abs_code_addr_immediate(int rdreg, int sizeshift, W64 imm) {
  put(TransOp(OP_add, rdreg, REG_trace, REG_imm, REG_zero, sizeshift, signext64(imm, 48)));
}

int TraceDecoder::bias_by_segreg(int basereg) {
  if (prefixes & (PFX_CS | PFX_DS | PFX_ES | PFX_FS | PFX_GS | PFX_SS)) {
    const SegmentRegister segment = (prefixes & PFX_FS)   ? SegmentRegister::fs
                                    : (prefixes & PFX_GS) ? SegmentRegister::gs
                                    : (prefixes & PFX_DS) ? SegmentRegister::ds
                                    : (prefixes & PFX_SS) ? SegmentRegister::ss
                                    : (prefixes & PFX_ES) ? SegmentRegister::es
                                                          : SegmentRegister::cs;

    int varoffs = offsetof_(Context, segment_bases) + segment_register_index(segment) * sizeof(address_t);

    TransOp ldp(OP_ld, REG_temp6, REG_ctx, REG_imm, REG_zero, 3, varoffs);
    ldp.internal = 1;
    put(ldp);
    put(TransOp(OP_add, REG_temp6, REG_temp6, basereg, REG_zero, 3));
    return REG_temp6;
  }

  return basereg;
}

void TraceDecoder::address_generate_and_load_or_store(int destreg, int srcreg, const DecodedOperand& memref, int opcode,
                                                      int datatype, int cachelevel, bool force_seg_bias, bool rmw) {
  //
  // In the address generation form used by internally generated
  // uops, we need the full virtual address, including the segment base
  // bias. However, technically LEA is not supposed to add the segment
  // base, even if prefixes implying a base are applied.
  //
  // Therefore, any internal uses that need the full virtual address
  // can use the special hint force_seg_bias to generate it,
  // but not actually include any load or store uops.
  //
  bool memop = isload(opcode) | isstore(opcode) | (opcode == OP_ld_pre);
  force_seg_bias |= memop;

  int imm_bits = (memop) ? 32 : 64;

  int basereg = arch_pseudo_reg_to_arch_reg[memref.mem.basereg];
  int indexreg = arch_pseudo_reg_to_arch_reg[memref.mem.indexreg];
  // ld rd = ra,rb,rc
  // ra = base
  // rb = offset or imm8
  // rc = reg to merge low bytes with

  //
  // Encoding rules for loads:
  //
  // Loads have only one addressing mode: basereg + immN, where
  // N is a small number of bits. The immediate is shifted left
  // by ld.size to allow for a greater range, assuming the immN
  // is always a multiple of the load data size. If immN is not
  // properly aligned, or it exceeds the field width, it cannot
  // be represented in the load and must be encoded outside as
  // a separate move immediate uop.
  //
  // Note that the rbimm field in the load uop is not actually
  // modified; this is for simulation purposes only, since real
  // microprocessors do not have unlimited immediate lengths.
  //

  bool imm_is_not_encodable = (lowbits(memref.mem.offset, memref.mem.size) != 0) |
                              (!fits_in_signed_nbit(memref.mem.offset >> memref.mem.size, imm_bits));

  if unlikely (opcode == OP_add) {
    // LEA and the like are always encodable since it's just an ADD uop:
    imm_is_not_encodable = 0;
  }

  W64s offset = memref.mem.offset;

  bool locked = ((prefixes & PFX_LOCK) != 0);

  if (basereg == REG_rip) {
    // [rip + imm32]: index always is zero and scale is 1
    // This mode is only possible in x86-64 code
    basereg = REG_zero;
    if (force_seg_bias)
      basereg = bias_by_segreg(basereg);

    if (memop) {
      abs_code_addr_immediate(REG_temp8, 3, Waddr(rip) + offset);
      TransOp ldst(opcode, destreg, REG_temp8, REG_imm, srcreg, memref.mem.size, 0);
      ldst.datatype = datatype;
      ldst.cachelevel = cachelevel;
      ldst.locked = locked;
      ldst.extshift = 0;
      put(ldst);
    } else {
      assert(opcode == OP_add);
      abs_code_addr_immediate(destreg, 3, Waddr(rip) + offset);
    }
  } else if (indexreg == REG_zero) {
    // [ra + imm32] or [ra]
    if (force_seg_bias)
      basereg = bias_by_segreg(basereg);
    if (imm_is_not_encodable) {
      put(TransOp(OP_add, REG_temp8, basereg, REG_imm, REG_zero, 3, offset));
      basereg = REG_temp8;
      offset = 0;
    }

    TransOp ldst(opcode, destreg, basereg, REG_imm, srcreg, memref.mem.size, offset);
    if (memop) {
      ldst.datatype = datatype;
      ldst.cachelevel = cachelevel;
      ldst.locked = locked;
      ldst.extshift = 0; // rmw;
    }
    put(ldst);
  } else if (offset == 0) {
    // [ra + rb*scale] or [rb*scale]
    if (force_seg_bias)
      basereg = bias_by_segreg(basereg);

    int tempreg = (memop) ? REG_temp8 : destreg;

    if (memref.mem.scale) {
      TransOp addop(OP_adda, tempreg, basereg, REG_zero, indexreg, (memop) ? 3 : memref.mem.size);
      addop.extshift = memref.mem.scale;
      put(addop);
    } else {
      put(TransOp(OP_add, tempreg, basereg, indexreg, REG_zero, (memop) ? 3 : memref.mem.size));
    }

    TransOp ldst(opcode, destreg, tempreg, REG_imm, srcreg, memref.mem.size, 0);
    if (memop) {
      // No need for this when we're only doing address generation:
      ldst.datatype = datatype;
      ldst.cachelevel = cachelevel;
      ldst.locked = locked;
      ldst.extshift = 0; // rmw;
    }
    put(ldst);
  } else {
    // [ra + imm32 + rb*scale]
    if (force_seg_bias)
      basereg = bias_by_segreg(basereg);

    if (imm_is_not_encodable) {
      put(TransOp(OP_add, REG_temp8, basereg, REG_imm, REG_zero, 3, offset));
      basereg = REG_temp8;
      offset = 0;
    }

    TransOp addop(OP_adda, REG_temp8, basereg, REG_zero, indexreg, 3);
    addop.extshift = memref.mem.scale;
    put(addop);
    TransOp ldst(opcode, destreg, REG_temp8, REG_imm, srcreg, memref.mem.size, offset);
    if (memop) {
      ldst.datatype = datatype;
      ldst.cachelevel = cachelevel;
      ldst.locked = locked;
      ldst.extshift = 0; // rmw;
    }
    put(ldst);
  }
}

void TraceDecoder::operand_load(int destreg, const DecodedOperand& memref, int opcode, int datatype, int cachelevel,
                                bool rmw) {
  address_generate_and_load_or_store(destreg, REG_zero, memref, opcode, datatype, cachelevel, false, rmw);
}

void TraceDecoder::result_store(int srcreg, int tempreg, const DecodedOperand& memref, int opcode, int datatype,
                                bool rmw) {
  address_generate_and_load_or_store(REG_mem, srcreg, memref, opcode, datatype, 0, 0, rmw);
}

void TraceDecoder::alu_reg_or_mem(int opcode, const DecodedOperand& rd, const DecodedOperand& ra, W32 setflags,
                                  int rcreg, bool flagsonly, bool isnegop, bool ra_rb_imm_form,
                                  W64s ra_rb_imm_form_rbimm) {
  if (flagsonly)
    prefixes &= ~PFX_LOCK;

  if ((rd.type == OPTYPE_REG) && ((ra.type == OPTYPE_REG) || (ra.type == OPTYPE_IMM))) {
    //
    // reg,reg or reg,imm
    //
    prefixes &= ~PFX_LOCK; // No locking on reg,reg
    assert(rd.reg.reg >= 0 && rd.reg.reg < APR_COUNT);
    if (ra.type == OPTYPE_REG)
      assert(ra.reg.reg >= 0 && ra.reg.reg < APR_COUNT);
    bool isimm = (ra.type == OPTYPE_IMM);
    int destreg = arch_pseudo_reg_to_arch_reg[rd.reg.reg];
    int srcreg = (isimm) ? REG_imm : arch_pseudo_reg_to_arch_reg[ra.reg.reg];
    int sizeshift = reginfo[rd.reg.reg].sizeshift;

    bool rdhigh = reginfo[rd.reg.reg].hibyte;
    bool rahigh = (isimm) ? 0 : reginfo[ra.reg.reg].hibyte;

    int rareg = destreg;
    if (rdhigh) {
      put(TransOp(OP_maskb, REG_temp2, REG_zero, rareg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      rareg = REG_temp2;
    }

    int rbreg = srcreg;
    if (rahigh) {
      put(TransOp(OP_maskb, REG_temp3, REG_zero, srcreg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      rbreg = REG_temp3;
    }

    //
    // Special case to break dependency chain for common idiom "xor X,X" => "xor zero,zero" => zero (always zero result when size >= 4 bytes)
    //
    if unlikely ((opcode == OP_xor) && (destreg == srcreg) && (sizeshift >= 2)) {
      assert(!(rahigh | rdhigh));
      rareg = REG_zero;
      rbreg = REG_zero;
    }

    if (flagsonly) {
      put(TransOp(opcode, REG_temp0, rareg, rbreg, rcreg, sizeshift, (isimm) ? ra.imm.imm : 0, 0, setflags));
    } else {
      if (isnegop) {
        rbreg = rareg;
        rareg = REG_zero;
      }
      if (ra_rb_imm_form) {
        put(TransOp(opcode, destreg, srcreg, REG_imm, (sizeshift >= 2) ? REG_zero : destreg, sizeshift,
                    ra_rb_imm_form_rbimm, 0, setflags));
      } else {
        if unlikely (isnegop && (sizeshift <= 1)) {
          // In decode-fast, neg r1 optimized to be => r1 = 0 - rd
          // but that only works for 32 and 64 bit operands
          // For 8 and 16 bit operands, use:
          // t0 = 0 - rd
          // mov rd=rd,t0
          put(TransOp(OP_sub, REG_temp0, REG_zero, destreg, REG_zero, sizeshift, 0, 0, setflags));
          put(TransOp(OP_mov, destreg, destreg, REG_temp0, REG_zero, sizeshift));
        } else {
          put(TransOp(opcode, (rdhigh) ? REG_temp2 : destreg, rareg, rbreg, rcreg, sizeshift, (isimm) ? ra.imm.imm : 0,
                      0, setflags));
          if (rdhigh) {
            put(TransOp(OP_maskb, destreg, destreg, REG_temp2, REG_imm, 3, 0, MaskControlInfo(56, 8, 56)));
          }
        }
      }
    }
  } else if ((rd.type == OPTYPE_REG) && (ra.type == OPTYPE_MEM)) {
    assert(rd.reg.reg >= 0 && rd.reg.reg < APR_COUNT);
    assert(ra.mem.basereg >= 0 && ra.mem.basereg < APR_COUNT);
    assert(ra.mem.indexreg >= 0 && ra.mem.indexreg < APR_COUNT);
    assert(ra.mem.scale >= 0 && ra.mem.scale <= 3);

    //
    // reg,[mem]
    //
    prefixes &= ~PFX_LOCK; // No locking on reg,[mem]
    int destreg = arch_pseudo_reg_to_arch_reg[rd.reg.reg];
    operand_load(REG_temp0, ra);

    bool rdhigh = reginfo[rd.reg.reg].hibyte;

    int rareg = destreg;
    if (rdhigh) {
      put(TransOp(OP_maskb, REG_temp2, REG_zero, destreg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      rareg = REG_temp2;
    }

    int sizeshift = reginfo[rd.reg.reg].sizeshift;
    if (flagsonly) {
      put(TransOp(opcode, REG_temp0, rareg, REG_temp0, rcreg, sizeshift, 0, 0, setflags));
    } else {
      if (ra_rb_imm_form) {
        put(TransOp(opcode, destreg, REG_temp0, REG_imm, (sizeshift >= 2) ? REG_zero : destreg, sizeshift,
                    ra_rb_imm_form_rbimm, 0, setflags));
      } else {
        put(TransOp(opcode, (rdhigh) ? REG_temp2 : destreg, rareg, REG_temp0, rcreg, sizeshift, 0, 0, setflags));
        if (rdhigh)
          put(TransOp(OP_maskb, destreg, destreg, REG_temp2, REG_imm, 3, 0, MaskControlInfo(56, 8, 56)));
      }
    }
  } else if ((rd.type == OPTYPE_MEM) && ((ra.type == OPTYPE_REG) || (ra.type == OPTYPE_IMM))) {
    //
    // [mem],reg or [mem],imm (rmw)
    //
    assert(rd.mem.basereg >= 0 && rd.mem.basereg < APR_COUNT);
    assert(rd.mem.indexreg >= 0 && rd.mem.indexreg < APR_COUNT);
    assert(rd.mem.scale >= 0 && rd.mem.scale <= 3);
    if (ra.type == OPTYPE_REG)
      assert(ra.reg.reg >= 0 && ra.reg.reg < APR_COUNT);

    bool isimm = (ra.type == OPTYPE_IMM);
    int srcreg = (isimm) ? REG_imm : arch_pseudo_reg_to_arch_reg[ra.reg.reg];
    operand_load(REG_temp0, rd, OP_ld, 0, 0, 1);

    int sizeshift = rd.mem.size;
    bool rahigh = (isimm) ? 0 : reginfo[ra.reg.reg].hibyte;

    if (rahigh) {
      put(TransOp(OP_maskb, REG_temp2, REG_zero, srcreg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      srcreg = REG_temp2;
    }

    if (isimm) {
      put(TransOp(opcode, REG_temp0, REG_temp0, REG_imm, rcreg, sizeshift, ra.imm.imm, 0, setflags));
      if (!flagsonly)
        result_store(REG_temp0, REG_temp3, rd, OP_st, 0, 1);
    } else {
      put(TransOp(opcode, REG_temp0, REG_temp0, srcreg, rcreg, sizeshift, 0, 0, setflags));
      if (!flagsonly)
        result_store(REG_temp0, REG_temp3, rd, OP_st, 0, 1);
    }
  } else if ((rd.type == OPTYPE_MEM) && (ra.type == OPTYPE_MEM)) {
    //
    // unary operations only: [mem],[samemem]
    //
    assert(rd.mem.basereg >= 0 && ra.mem.basereg < APR_COUNT);
    assert(rd.mem.indexreg >= 0 && ra.mem.indexreg < APR_COUNT);
    assert(rd.mem.scale >= 0 && ra.mem.scale <= 3);

    operand_load(REG_temp0, rd);
    int sizeshift = rd.mem.size;
    put(TransOp(opcode, REG_temp0, (isnegop) ? REG_zero : REG_temp0, REG_temp0, rcreg, sizeshift, 0, 0, setflags));
    if (!flagsonly)
      result_store(REG_temp0, REG_temp3, rd);
  }

  if unlikely (no_partial_flag_updates_per_insn && (setflags != (SETFLAG_ZF | SETFLAG_CF | SETFLAG_OF))) {
    put(TransOp(OP_collcc, REG_temp10, REG_zf, REG_cf, REG_of, 3, 0, 0, FLAGS_DEFAULT_ALU));
  }
}

void TraceDecoder::move_reg_or_mem(const DecodedOperand& rd, const DecodedOperand& ra, int force_rd) {
  prefixes &= ~PFX_LOCK; // No locking unless both src and dest are memory

  if ((rd.type == OPTYPE_REG) && ((ra.type == OPTYPE_REG) || (ra.type == OPTYPE_IMM))) {
    //
    // reg,reg
    //

    bool isimm = (ra.type == OPTYPE_IMM);
    int destreg = (force_rd != REG_zero) ? force_rd : arch_pseudo_reg_to_arch_reg[rd.reg.reg];
    int srcreg = (isimm) ? REG_imm : arch_pseudo_reg_to_arch_reg[ra.reg.reg];
    int sizeshift = reginfo[rd.reg.reg].sizeshift;

    bool rdhigh = (force_rd != REG_zero) ? false : reginfo[rd.reg.reg].hibyte;
    bool rahigh = (isimm) ? 0 : reginfo[ra.reg.reg].hibyte;

    if (rdhigh || rahigh) {
      int maskctl = (rdhigh && !rahigh) ? MaskControlInfo(56, 8, 56) : // insert high byte
                        (!rdhigh && rahigh) ? MaskControlInfo(0, 8, 8)
                                            : // extract high byte
                        (rdhigh && rahigh) ? MaskControlInfo(56, 8, 0)
                                           :      // move between high bytes
                        MaskControlInfo(0, 8, 0); // move between low bytes
      put(TransOp(OP_maskb, destreg, destreg, srcreg, REG_imm, 3, (isimm) ? ra.imm.imm : 0, maskctl));
    } else {
      // must be at least 16 bits
      // On x86-64, only 8-bit and 16-bit ops need to be merged; 32-bit is zero extended to full 64 bits:
      put(TransOp(OP_mov, destreg, (sizeshift < 2) ? destreg : REG_zero, srcreg, REG_zero, sizeshift,
                  (isimm) ? ra.imm.imm : 0));
    }
  } else if ((rd.type == OPTYPE_REG) && (ra.type == OPTYPE_MEM)) {
    //
    // reg,[mem]
    //
    int destreg = (force_rd != REG_zero) ? force_rd : arch_pseudo_reg_to_arch_reg[rd.reg.reg];
    int sizeshift = reginfo[rd.reg.reg].sizeshift;

    if ((sizeshift >= 2) || (force_rd != REG_zero)) {
      // zero extend 32-bit to 64-bit or just load as 64-bit:
      operand_load(destreg, ra);
    } else {
      // need to merge 8-bit or 16-bit data:
      operand_load(REG_temp0, ra);
      if (reginfo[rd.reg.reg].hibyte)
        put(TransOp(OP_maskb, destreg, destreg, REG_temp0, REG_imm, 3, 0, MaskControlInfo(56, 8, 56)));
      else
        put(TransOp(OP_mov, destreg, destreg, REG_temp0, REG_zero, sizeshift));
    }
  } else if ((rd.type == OPTYPE_MEM) && ((ra.type == OPTYPE_REG) || (ra.type == OPTYPE_IMM))) {
    //
    // [mem],reg
    //
    bool isimm = (ra.type == OPTYPE_IMM);
    int srcreg = (isimm) ? REG_imm : arch_pseudo_reg_to_arch_reg[ra.reg.reg];

    bool rahigh = (isimm) ? 0 : reginfo[ra.reg.reg].hibyte;
    if (isimm) {
      // We need to load the immediate separately in any case since stores do not accept immediates:
      put(TransOp(OP_mov, REG_temp1, REG_zero, REG_imm, REG_zero, 3, ra.imm.imm));
      result_store(REG_temp1, REG_temp0, rd);
    } else if (rahigh) {
      put(TransOp(OP_maskb, REG_temp1, REG_zero, srcreg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      result_store(REG_temp1, REG_temp0, rd);
    } else {
      result_store(srcreg, REG_temp0, rd);
    }
  }
}

void TraceDecoder::signext_reg_or_mem(const DecodedOperand& rd, DecodedOperand& ra, int rasize, bool zeroext) {
  int rdreg = arch_pseudo_reg_to_arch_reg[rd.reg.reg];
  int rdsize = reginfo[rd.reg.reg].sizeshift;
  prefixes &= ~PFX_LOCK; // No locking unless both src and dest are memory

  if ((rd.type == OPTYPE_REG) && (ra.type == OPTYPE_REG)) {
    //
    // reg,reg
    //

    int rareg = arch_pseudo_reg_to_arch_reg[ra.reg.reg];

    assert(!reginfo[rd.reg.reg].hibyte);
    bool rahigh = reginfo[ra.reg.reg].hibyte;

    // For movsx, it is not possible to have rd be 8 bits: it's always 16/32/64 bits.
    if (rahigh) {
      put(TransOp(OP_maskb, REG_temp0, REG_zero, rareg, REG_imm, 3, 0, MaskControlInfo(0, 8, 8)));
      rareg = REG_temp0;
    }

    assert(rasize < 3); // must be at most 32 bits
    // On x86-64, only 8-bit and 16-bit ops need to be merged; 32-bit is zero extended to full 64 bits:
    if (zeroext && rdsize >= 2) {
      // Just use regular move
      put(TransOp(OP_mov, rdreg, REG_zero, rareg, REG_zero, rasize));
    } else {
      TransOp transop(OP_maskb, rdreg, (rdsize < 2) ? rdreg : REG_zero, rareg, REG_imm, rdsize, 0,
                      MaskControlInfo(0, (1 << rasize) * 8, 0));
      transop.cond = (zeroext) ? 1 : 2;
      put(transop);
    }
  } else if ((rd.type == OPTYPE_REG) && (ra.type == OPTYPE_MEM)) {
    //
    // reg,[mem]
    //
    ra.mem.size = rasize;

    if (rdsize >= 2) {
      bool signext_to_32bit = (rdsize == 2) & (!zeroext);
      // zero extend 32-bit to 64-bit or just load as 64-bit:
      operand_load((signext_to_32bit) ? REG_temp8 : rdreg, ra, (zeroext) ? OP_ld : OP_ldx);
      // sign extend and then zero high 32 bits (old way was ldxz uop):
      if (signext_to_32bit)
        put(TransOp(OP_mov, rdreg, REG_zero, REG_temp8, REG_zero, 2));
    } else {
      // need to merge 8-bit or 16-bit data:
      operand_load(REG_temp0, ra, (zeroext) ? OP_ld : OP_ldx);
      put(TransOp(OP_mov, rdreg, rdreg, REG_temp0, REG_zero, rdsize));
    }
  }
}

void TraceDecoder::microcode_assist(int assistid, Waddr selfrip, Waddr nextrip) {
  used_microcode_assist = 1;
  abs_code_addr_immediate(REG_selfrip, 3, (Waddr)selfrip);
  abs_code_addr_immediate(REG_nextrip, 3, (Waddr)nextrip);
  if (!last_flags_update_was_atomic)
    put(TransOp(OP_collcc, REG_temp0, REG_zf, REG_cf, REG_of, 3, 0, 0, FLAGS_DEFAULT_ALU));
  TransOp transop(OP_brp, REG_rip, REG_zero, REG_zero, REG_zero, 3);
  transop.riptaken = transop.ripseq = (Waddr)assistid;
  bb.rip_taken = assistid;
  bb.rip_not_taken = assistid;
  bb.brtype = BRTYPE_BARRIER;
  put(transop);
}

//
// Core Decoder
//
void TraceDecoder::decode_prefixes() {
  prefixes = 0;
  rex = 0;

  for (;;) {
    byte b = insnbytes[byteoffset];
    W32 prefix = (use64) ? prefix_map_x86_64[b] : prefix_map_x86[b];
    if (!prefix)
      break;
    if (rex) {
      // REX is ignored when followed by another prefix:
      rex = 0;
      prefixes &= ~PFX_REX;
    }
    prefixes |= prefix;
    if (prefix == PFX_REX) {
      rex = b;
    }
    byteoffset++;
    rip++;
  }
}

void TraceDecoder::split(bool after) {
  Waddr target = (after) ? rip : ripstart;
  if (!after)
    assert(!first_insn_in_bb());

  //
  // Append to the previous insn (no pre-flush) if split before
  // the current insn. This forces the total x86 insn count to
  // match the real code (i.e. do not count branch caps):
  //
  put(TransOp(OP_collcc, REG_temp0, REG_zf, REG_cf, REG_of, 3, 0, 0, FLAGS_DEFAULT_ALU));

  TransOp br(OP_bru, REG_rip, REG_zero, REG_zero, REG_zero, 3);
  br.riptaken = target;
  br.ripseq = target;
  put(br);
  bb.rip_taken = target;
  bb.rip_not_taken = target;
  bb.brtype = BRTYPE_SPLIT;
  join_with_prev_insn = (!after);
  end_of_block = 1;
}

void print_invalid_insns(MachineImpl& machine, int op, const byte* ripstart, const byte* rip, int valid_byte_count,
                         const PageFaultErrorCode& pfec, Waddr faultaddr) {
  if (pfec) {
    logging::println(
        logging::INFO,
        "translate: page fault at iteration {}, {}, commits: ripstart {}, rip {}: required {} more bytes but "
        "only fetched {} bytes; page fault error code: {}",
        machine.iterations, machine.total_user_insns_committed, (const void*)ripstart, (const void*)rip,
        (rip - ripstart), valid_byte_count, pfec);
    logging::flush();
  } else {
    logging::println(
        logging::INFO,
        "translate: invalid opcode at iteration {}: {} commits (at ripstart {}, rip {}); may be speculative",
        machine.iterations, machine.total_user_insns_committed, (const void*)ripstart, (const void*)rip);
    logging::flush();
#if 0
    if (!config.dumpcode_filename.empty()) {
      byte insnbuf[256];
      PageFaultErrorCode copypfec;
      int valid_byte_count = contextof(0).copy_from_user(insnbuf, (Waddr)rip, sizeof(insnbuf), copypfec, faultaddr);
      odstream os(config.dumpcode_filename);
      os.write(insnbuf, sizeof(insnbuf));
      os.close();
    }
#endif
  }
}

void assist_invalid_opcode(Context& ctx) {
  ctx.commitarf[REG_rip] = ctx.commitarf[REG_selfrip];
  ctx.propagate_x86_exception(EXCEPTION_x86_invalid_opcode);
}

static const bool log_code_page_ops = 0;

BasicBlockCache::~BasicBlockCache() {
  flush();
}

bool BasicBlockCache::invalidate(BasicBlock* bb, int reason) {
  BasicBlockChunkList* pagelist;
  if unlikely (bb->refcount) {
    logging::println(logging::INFO, "Warning: basic block {} <bb> is still in use somewhere (refcount {})", (void*)bb,
                     bb->refcount);
    return false;
  }

  if unlikely (bbcache_dump_file) {
    if (std::fwrite(static_cast<BasicBlockBase*>(bb), sizeof(BasicBlockBase), 1, bbcache_dump_file) != 1 ||
        std::fwrite(bb->transops, sizeof(TransOp), bb->count, bbcache_dump_file) != (size_t)bb->count) {
      logging::println(logging::WARNING, "Failed to write basic block {} to bb dump file", (void*)bb);
    }
  }

  pagelist = pages_.get(bb->rip.mfnlo);
  logging::println(logging::INFO, "Remove bb {} ({}, {} bytes) from low page list {}: loc {}:{}", (void*)bb,
                   (void*)(Waddr)bb->rip, bb->bytes, (void*)pagelist, (void*)bb->mfnlo_loc.chunk, bb->mfnlo_loc.index);
  assert(pagelist);
  pagelist->remove(bb->mfnlo_loc);

  int page_crossing = ((lowbits(bb->rip, 12) + (bb->bytes - 1)) >> 12);
  if (page_crossing) {
    pagelist = pages_.get(bb->rip.mfnhi);
    logging::println(logging::DEBUG, "Remove bb {} ({}, {} bytes) from high page list {}: loc {}:{}", (void*)bb,
                     (void*)(Waddr)bb->rip, bb->bytes, (void*)pagelist, (void*)bb->mfnhi_loc.chunk,
                     bb->mfnhi_loc.index);
    assert(pagelist);
    pagelist->remove(bb->mfnhi_loc);
  }

  remove(bb);
  stats.decoder.bbcache.count = this->count;
  stats.decoder.bbcache.invalidates[reason]++;

  bb->free();
  return true;
}

bool BasicBlockCache::invalidate(const RIPVirtPhys& rvp, int reason) {
  BasicBlock* bb = get(rvp);
  if (!bb)
    return true;
  return invalidate(bb, reason);
}

//
// Find the number of cached BBs on a given physical page
//
int BasicBlockCache::get_page_bb_count(Waddr mfn) {
  if unlikely (mfn == RIPVirtPhys::INVALID)
    return 0;

  BasicBlockChunkList* pagelist = pages_.get(mfn);

  if unlikely (!pagelist)
    return 0;

  return pagelist->count();
}

//
// Invalidate any basic blocks on the specified physical page.
// This function is suitable for calling from a reclaim handler
// when we run out of memory (it may will allocate any memory).
//
bool BasicBlockCache::invalidate_page(AddressSpace& asp, Waddr mfn, int reason) {
  //
  // We may try to invalidate the special invalid mfn if SMC
  // occurs on a page where the high virtual page is invalid.
  //
  if unlikely (mfn == RIPVirtPhys::INVALID)
    return 0;

  BasicBlockChunkList* pagelist = pages_.get(mfn);

  if (log_code_page_ops)
    logging::println(logging::INFO, "Invalidate page mfn {}: pagelist {} has {} entries (dirty? {})", mfn,
                     (void*)pagelist, (pagelist ? pagelist->count() : 0), asp.isdirty(mfn));

  asp.cleardirty(mfn);

  /* In the event that there are no basicblocks for this page,
   * i.e. we have no pagelist, then this is a successfull NOP
   */
  if unlikely (!pagelist) {
    return true;
  }

  /* invalidate calls pagelist.remove so we can't really use an iterator
   * Hack around this by resetting the iterator after every call to invalidate
   * Reusing the same iterator should atleast keep the performance impact small
   */
  BasicBlockChunkList::Iterator iter(pagelist);
  while (pagelist->count() > 0) {
    iter.reset(pagelist);

    BasicBlock* bb = *iter.next();
    if (log_code_page_ops) {
      logging::println(logging::DEBUG, "  Invalidate bb {} ({}, {} bytes)", (void*)bb, (void*)(Waddr)bb->rip,
                       bb->bytes);
    }

    if unlikely (!invalidate(bb, reason)) {
      if (log_code_page_ops) {
        logging::println(logging::DEBUG, "  Could not invalidate bb {} ({}, {} bytes): still has refcount {}",
                         (void*)bb, bb->rip, bb->bytes, bb->refcount);
      }
      return false;
    }
  }

  /* assert that a call to pagelist->clear() is unecessary */
  assert(pagelist->head == NULL);
  assert(pagelist->elemcount == 0);

  stats.decoder.pagecache.count = pages_.count;
  stats.decoder.pagecache.invalidates[reason]++;

  return true;
}

//
// Scan through the BB cache and try to free up
// <bytesreq> bytes, starting with the least
// recently used BBs.
//
int BasicBlockCache::reclaim(size_t bytesreq, int urgency) {
  if (!count)
    return 0;

  logging::println(logging::DEBUG, "Reclaiming cached basic blocks:");

  stats.decoder.reclaim_rounds++;

  W64 oldest = std::numeric_limits<W64>::max();
  W64 newest = 0;
  W64 average = 0;

  int n = 0;

  Iterator iter(this);
  BasicBlock* bb;

  while ((bb = iter.next())) {
    oldest = std::min(oldest, bb->lastused);
    newest = std::max(newest, bb->lastused);
    average += bb->lastused;
    n++;
  }

  assert(count == n);
  assert(n > 0);
  average /= n;

  if unlikely (urgency >= /*MAX_URGENCY*/ 10000000) {
    //
    // The allocator is so strapped for memory, we need to free
    // everything possible at all costs:
    //
    average = std::numeric_limits<W64s>::max();
  }

  logging::println(logging::DEBUG, "Before:");
  logging::println(logging::DEBUG, "  Basic blocks:   {:>12}", count);
  logging::println(logging::DEBUG, "  Oldest cycle:   {:>12}", oldest);
  logging::println(logging::DEBUG, "  Average cycle:  {:>12}", average);
  logging::println(logging::DEBUG, "  Newest cycle:   {:>12}", newest);

  //
  // Reclaim all objects older than the average
  //
  n = 0;
  int reclaimed_objs = 0;

  iter.reset(this);
  while ((bb = iter.next())) {
    if unlikely (bb->refcount) {
      //
      // We cannot invalidate anything that's still in the pipeline.
      // If this is required, the pipeline must be flushed before
      // the forced invalidation can occur.
      //
      logging::println(logging::INFO, "Warning: eligible bb {} {} (lastused {}) still has refcount {}", (void*)bb,
                       (void*)(Waddr)bb->rip, bb->lastused, bb->refcount);
      continue;
    }

    // We use '<=' to guarantee even a uniform distribution will eventually be reclaimed:
    if likely (bb->lastused <= average) {
      reclaimed_objs++;
      invalidate(bb, INVALIDATE_REASON_RECLAIM);
    }
    n++;
  }

  logging::println(logging::DEBUG, "After:");
  logging::println(logging::DEBUG, "  Basic blocks:   {:>12} BBs reclaimed", reclaimed_objs);
  logging::println(logging::DEBUG, "  New pool size:  {:>12} BBs", count);
  logging::flush();


  //
  // Reclaim per-page chunklist heads
  //

  {
    BasicBlockPageCache::Iterator iter(&pages_);
    BasicBlockChunkList* page;
    int pages_freed = 0;

    logging::println(logging::DEBUG, "Scanning {} code pages:", pages_.count);
    while ((page = iter.next())) {
      if (page->empty()) {
        if (!page->refcount) {
          logging::println(logging::DEBUG, "  mfn {} has no entries; freeing", page->mfn);
          pages_.remove(page);
          delete page;
          pages_freed++;
        } else {
          logging::println(logging::DEBUG, "  mfn {} still has refs to it: cannot free it yet", page->mfn);
        }
      }
    }

    logging::println(logging::DEBUG, "Freed {} empty pages", pages_freed);
  }

  return n;
}

//
// Flush the entire basic block cache immediately.
// All basic blocks are flushed: no remaining
// references are allowed.
//
void BasicBlockCache::flush() {
  logging::println(logging::DEBUG, "Flushing basic block cache:");

  stats.decoder.reclaim_rounds++;

  {
    Iterator iter(this);
    BasicBlock* bb;
    while ((bb = iter.next())) {
      invalidate(bb, INVALIDATE_REASON_RECLAIM);
    }
  }

  //
  // Reclaim per-page chunklist heads
  //

  {
    BasicBlockPageCache::Iterator iter(&pages_);
    BasicBlockChunkList* page;
    int pages_freed = 0;

    while ((page = iter.next())) {
      assert(page->empty());
      pages_.remove(page);
      delete page;
    }
  }
}

void assist_exec_page_fault(Context& ctx) {
  //
  // We need to check if faultaddr is now a valid page, since the page tables
  // could have been updated since the cut-off basic block was speculatively
  // translated, such that the block is now valid. If the page at faultaddr
  // in REG_ar1 is now valid, we need to invalidate the currently executing
  // translation and tell the main loop to start translating again at the
  // cut-off instruction's starting byte.
  //
  Waddr faultaddr = ctx.commitarf[REG_ar1];
  PageFaultErrorCode pfec = ctx.commitarf[REG_ar2];

  bool page_now_valid = ctx.address_space->check(faultaddr, Protection::execute);
  if unlikely (page_now_valid) {
    logging::println(
        logging::WARNING,
        "Spurious PageFaultOnExec detected at fault rip {} with faultaddr {} @ {} user commits ({} cycles)",
        (void*)(Waddr)ctx.commitarf[REG_selfrip], (void*)faultaddr, ctx.machine_impl->total_user_insns_committed,
        ctx.machine_impl->sim_cycle);
    ctx.machine_impl->bbcache->invalidate(RIPVirtPhys(ctx.commitarf[REG_selfrip]).update(ctx),
                                          INVALIDATE_REASON_SPURIOUS);
    ctx.commitarf[REG_rip] = ctx.commitarf[REG_selfrip];
    return;
  }

  ctx.commitarf[REG_rip] = ctx.commitarf[REG_selfrip];
  ctx.propagate_x86_exception(EXCEPTION_x86_page_fault, pfec, faultaddr);
}

void assist_gp_fault(Context& ctx) {
  ctx.commitarf[REG_rip] = ctx.commitarf[REG_selfrip];
  ctx.propagate_x86_exception(EXCEPTION_x86_gp_fault, ctx.commitarf[REG_ar1]);
}

bool TraceDecoder::invalidate() {
  if likely ((rip - bb.rip) > valid_byte_count) {
    int mfn = 0;
    logging::println(
        logging::WARNING,
        "Translation crosses into invalid page (mfn {}): ripstart {}, rip {}, faultaddr {}; expected {} bytes "
        "but only got {} (next page {})",
        mfn, (void*)ripstart, (void*)rip, (void*)faultaddr, (rip - ripstart), valid_byte_count,
        (void*)(Waddr)ceil(ripstart, 4096));

    if likely (split_invalid_basic_blocks && (!first_insn_in_bb())) {
      //
      // Split BB before the first invalid instruction, such that the
      // instruction in question is the first one in the next BB.
      // Outside code is responsible for dispatching exceptions.
      //
      invalid = 0;
      rip = ripstart;
      split_before();
      return false;
    } else {
      outcome = (faultaddr == bb.rip.rip) ? DECODE_OUTCOME_ENTRY_PAGE_FAULT : DECODE_OUTCOME_OVERLAP_PAGE_FAULT;
      print_invalid_insns(*machine, op, (const byte*)ripstart, (const byte*)rip, valid_byte_count, pfec, faultaddr);
      abs_code_addr_immediate(REG_ar1, 3, faultaddr);
      immediate(REG_ar2, 3, pfec);
      microcode_assist(ASSIST_EXEC_PAGE_FAULT, ripstart, faultaddr);
    }
  } else {
    // The instruction-specific decoder may have already set the outcome type
    if (outcome == DECODE_OUTCOME_OK)
      outcome = DECODE_OUTCOME_INVALID_OPCODE;

    logging::println(logging::WARNING, "Invalid opcode at {}: split_invalid_basic_blocks {}, first_insn_in_bb? {}",
                     (void*)ripstart, split_invalid_basic_blocks, first_insn_in_bb());
    print_invalid_insns(*machine, op, (const byte*)ripstart, (const byte*)rip, valid_byte_count, 0, faultaddr);

    if likely (split_invalid_basic_blocks && (!first_insn_in_bb())) {
      //
      // Split BB before the first invalid instruction, such that the
      // instruction in question is the first one in the next BB.
      // Outside code is responsible for dispatching exceptions.
      //
      invalid = 0;
      bb.invalidblock = 0;
      outcome = DECODE_OUTCOME_OK;
      rip = ripstart;
      split_before();
      return false;
    } else {
      switch (outcome) {
      case DECODE_OUTCOME_INVALID_OPCODE:
        microcode_assist(ASSIST_INVALID_OPCODE, ripstart, rip);
        break;
      case DECODE_OUTCOME_GP_FAULT:
        pfec = 0;
        microcode_assist(ASSIST_GP_FAULT, ripstart, rip);
        break;
      default:
        logging::println(logging::ERROR, "Unexpected decoder outcome: {}", outcome);
        break;
      }
    }
  }
  end_of_block = 1;
  user_insn_count++;
  bb.invalidblock = 1;
  flush();
  return true;
}

//
// Generate a memory fence of the specified type.
//
bool TraceDecoder::memory_fence_if_locked(bool end_of_x86_insn, int type) {
  if likely (!(prefixes & PFX_LOCK))
    return false;

  if (split_basic_block_at_locks_and_fences) {
    if (end_of_x86_insn) {
      //
      // Final mf that terminates the insn (always in its own BB): terminate the BB here.
      //
      TransOp mf(OP_mf, REG_temp0, REG_zero, REG_zero, REG_zero, 0);
      mf.extshift = type;
      put(mf);
      bb.mfence = 1;
      split_after();
      return true;
    } else {
      //
      // First uop in x86 insn
      //
      if (!first_insn_in_bb()) {
        //
        // This is not the first x86 insn in the BB, but it is the first uop in the insn.
        // Split the existing BB at this point, so we get a sequence like this:
        //
        // insn1
        // insn2
        // (insn3):
        //   mf
        //   bru   insn3
        //
        split_before();
        return true;
      } else {
        //
        // This is the first insn in the BB and the first mf uop in that insn:
        // just emit the mf and continue with the other uops in the insn
        //
        TransOp mf(OP_mf, REG_temp0, REG_zero, REG_zero, REG_zero, 0);
        mf.extshift = type;
        put(mf);
        bb.mfence = 1;
        return false; // continue emitting other uops
      }
    }
  } else {
    //
    // Always emit fence intermingled with other insns
    //
    TransOp mf(OP_mf, REG_temp0, REG_zero, REG_zero, REG_zero, 0);
    mf.extshift = type;
    put(mf);
  }

  // never reached:
  return false;
}

//
// Fill the insnbytes buffer as much as possible,
// properly handling x86 semantics where the insn
// extends onto an invalid page. Return the number
// of valid bytes, if any.
//
// We limit BBs to at most insnbytes_bufsize; this ensures
// we can do a bulk copy of the potentially unaligned
// instruction bytes into insnbytes once at the start,
// rather than having to constantly do checks.
//

int TraceDecoder::fillbuf(Context& ctx, byte* insnbytes, int insnbytes_bufsize) {
  machine = ctx.machine_impl;
  assert(machine);
  this->insnbytes = insnbytes;
  this->insnbytes_bufsize = insnbytes_bufsize;
  byteoffset = 0;
  faultaddr = 0;
  pfec = 0;
  invalid = 0;
  valid_byte_count = ctx.copy_from_user(insnbytes, bb.rip, insnbytes_bufsize, pfec, faultaddr, true, ptelo, ptehi);
  return valid_byte_count;
}

//
// Decode and translate one x86 instruction
//
bool TraceDecoder::translate() {
  opsize_prefix = 0;
  addrsize_prefix = 0;
  first_uop_in_insn = 1;
  ripstart = rip;

  decode_prefixes();

  {
    std::string prefix_str;
    foreach (i, PFX_count) {
      if (prefixes & (1 << i))
        prefix_str += std::format(" {}", prefix_names[i]);
    }
    logging::println(logging::DEBUG, "prefixes = {}:{}", prefixes, prefix_str);
  }

  if (prefixes & PFX_ADDR)
    addrsize_prefix = 1;

  bool uses_sse = 0;
  op = fetch1();
  bool need_modrm = onebyte_has_modrm[op];
  if (op == 0x0f) {
    op = fetch1();
    need_modrm = twobyte_has_modrm[op];

    if (twobyte_uses_SSE_prefix[op]) {
      uses_sse = 1;
      if (prefixes & PFX_DATA) // prefix byte 0x66, typically OPpd
        op |= 0x500;
      else if (prefixes & PFX_REPNZ) // prefix byte 0xf2, typically OPsd
        op |= 0x400;
      else if (prefixes & PFX_REPZ) // prefix byte 0xf3, typically OPss
        op |= 0x200;
      else
        op |= 0x300; // no prefix byte, typically OPps
    } else {
      op |= 0x100;
    }
  }

  // SSE uses 0x66 prefix for an opcode extension:
  if (!uses_sse && (prefixes & PFX_DATA))
    opsize_prefix = 1;

  modrm = ModRMByte((need_modrm) ? fetch1() : 0);

  if (inrange(op, 0xd8, 0xdf)) {
    op = 0x600 | (lowbits(op, 3) << 4) | modrm.reg;
  }

  bool rc;

  logging::println(logging::DEBUG, "Decoding op {} (class {}) @ {}", hexstring(op, 12), (op >> 8), (void*)ripstart);

  is_x87 = 0;
  is_sse = 0;

  invalid |= ((rip - (Waddr)bb.rip) > valid_byte_count);

  if (invalid) {
    invalidate();
    user_insn_count++;
    end_of_block = 1;
    flush();
    return false;
  }

  switch (op >> 8) {
  case 0:
  case 1: {
    rc = decode_fast();

    // Try again with the complex decoder if needed
    bool iscomplex = ((rc == 0) & (!invalid));
    some_insns_complex |= iscomplex;
    if (iscomplex)
      rc = decode_complex();

    if unlikely (used_microcode_assist) {
      stats.decoder.x86_decode_type[DECODE_TYPE_ASSIST]++;
    } else {
      stats.decoder.x86_decode_type[DECODE_TYPE_FAST] += (!iscomplex);
      stats.decoder.x86_decode_type[DECODE_TYPE_COMPLEX] += iscomplex;
    }

    break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
    stats.decoder.x86_decode_type[DECODE_TYPE_SSE]++;
    rc = decode_sse();
    break;
  case 6:
    stats.decoder.x86_decode_type[DECODE_TYPE_X87]++;
    rc = decode_x87();
    break;
  default: {
    assert(false);
  }
  } // switch

  if (!rc)
    return rc;

  user_insn_count++;

  assert(!invalid);

  if (end_of_block) {
    // Block ended with a branch: close the uop and exit
    stats.decoder.bb_decode_type.all_insns_fast += (!some_insns_complex);
    stats.decoder.bb_decode_type.some_complex_insns += some_insns_complex;
    flush();
    return false;
  } else {
    // Block did not end with a branch: do we have more room for another x86 insn?
    if ( // ((MAX_BB_UOPS - bb.count) < (MAX_TRANSOPS_PER_USER_INSN-2)) ||
        ((rip - bb.rip) >= (insnbytes_bufsize - 15)) || ((rip - bb.rip) >= valid_byte_count) ||
        (user_insn_count >= MAX_BB_X86_INSNS) || (rip == stop_at_rip)) {
      logging::println(logging::DEBUG, "Basic block {} too long: cutting at {} transops ({} currently in buffer)",
                       (void*)(Waddr)bb.rip, bb.count, transbufcount);
      // bb.rip_taken and bb.rip_not_taken were already filled out for the last instruction.
      if unlikely (!last_flags_update_was_atomic) {
        logging::println(logging::DEBUG, "Basic block {} had non-atomic flags update: adding collcc",
                         (void*)(Waddr)bb.rip);
        put(TransOp(OP_collcc, REG_temp0, REG_zf, REG_cf, REG_of, 3, 0, 0, FLAGS_DEFAULT_ALU));
      }
      split_after();
      stats.decoder.bb_decode_type.all_insns_fast += (!some_insns_complex);
      stats.decoder.bb_decode_type.some_complex_insns += some_insns_complex;
      flush();
      return false;
    } else {
      bool ok_to_continue = flush();
      return ok_to_continue;
    }
  }
#undef DECODE
}

std::string format_flags(W64 flags) {
  std::string result = std::format("0x{:08x} = [", flags);

  for (int i = (FLAG_COUNT - 1); i >= 0; i--) {
    if (bit(flags, i))
      result += std::format(" {}", x86_flag_names[i]);
    else
      result += " -";
  }
  result += " ] ";
  return result;
}

//
// Translate one basic block. This function always returns
// a BasicBlock, except in the very rare case where one or
// both covered mfns are dirty and must be invalidated, and
// the invalidation fails because some other object has
// references to some of the basic blocks.
//
BasicBlock* BasicBlockCache::translate(Context& ctx, const RIPVirtPhys& rvp) {
  AddressSpace& asp = *ctx.address_space;
  Options& config = ctx.machine_impl->config;

  if unlikely ((rvp.rip == config.log.start_log_at_rip) && (rvp.rip != 0xffffffffffffffffULL)) {
    config.log.start_log_at_iteration = 0;
    logenable = 1;
  }

  if unlikely (asp.isdirty(rvp.mfnlo)) {
    if (log_code_page_ops)
      logging::println(logging::TRACE, "Pre-invalidate low mfn for {}", rvp);
    if unlikely (!invalidate_page(asp, rvp.mfnlo, INVALIDATE_REASON_DIRTY))
      return null;
  }

  if unlikely (asp.isdirty(rvp.mfnhi)) {
    if (log_code_page_ops)
      logging::println(logging::TRACE, "Pre-invalidate high mfn for {}", rvp);
    if unlikely (!invalidate_page(asp, rvp.mfnhi, INVALIDATE_REASON_DIRTY))
      return null;
  }

  BasicBlock* bb = get(rvp);
  if likely (bb)
    return bb;

  byte insnbuf[MAX_BB_BYTES];

  TraceDecoder trans(rvp);
  trans.fillbuf(ctx, insnbuf, sizeof(insnbuf));

  if (log_code_page_ops) {
    logging::println(logging::INFO, "Translating {} ({} bytes valid) at {} cycles, {} commits", rvp,
                     trans.valid_byte_count, ctx.machine_impl->sim_cycle, ctx.machine_impl->total_user_insns_committed);
  }

  if (rvp.mfnlo == RIPVirtPhys::INVALID) {
    assert(trans.valid_byte_count == 0);
  }

  for (;;) {
    logging::println(logging::DEBUG, "rip {}, relrip = {}", (void*)trans.rip, (void*)(trans.rip - trans.bb.rip));
    if (!trans.translate())
      break;
  }

  trans.bb.hitcount = 0;
  trans.bb.predcount = 0;
  bb = trans.bb.clone();
  //
  // Acquire a reference to the new basic block right away,
  // since we make allocations below that might reclaim it
  // out from under us.
  //
  bb->acquire();

  add(bb);
  stats.decoder.bbcache.count = this->count;
  stats.decoder.bbcache.inserts++;

  stats.decoder.throughput.basic_blocks++;

  BasicBlockChunkList* pagelist;

  asp.cleardirty(bb->rip.mfnlo);
  pagelist = pages_.get(bb->rip.mfnlo);
  if (!pagelist) {
    pagelist = new BasicBlockChunkList(bb->rip.mfnlo);
    pagelist->refcount++;
    pages_.add(pagelist);
    stats.decoder.pagecache.inserts++;
    stats.decoder.pagecache.count = pages_.count;
    pagelist->refcount--;
  }
  pagelist->refcount++;

  //
  // WARNING: in extreme low memory situations, if this needs memory
  // (to add a new chunk to the list), the list itself may get freed
  // since it technically does not have any MFNs on it yet. We need
  // to somehow lock it with a refcount to prevent this.
  //
  pagelist->add(bb, bb->mfnlo_loc);
  if (log_code_page_ops)
    logging::println(logging::TRACE, "Add bb {} ({}, {} bytes) to low page list {}: loc {}:{}", (void*)bb, bb->rip,
                     bb->bytes, (void*)pagelist, (void*)bb->mfnlo_loc.chunk, bb->mfnlo_loc.index);

  int page_crossing = ((lowbits(bb->rip, 12) + (bb->bytes - 1)) >> 12);

  if (page_crossing) {
    asp.cleardirty(bb->rip.mfnhi);
    BasicBlockChunkList* pagelisthi = pages_.get(bb->rip.mfnhi);
    if (!pagelisthi) {
      pagelisthi = new BasicBlockChunkList(bb->rip.mfnhi);
      pagelisthi->refcount++;
      pages_.add(pagelisthi);
      stats.decoder.pagecache.inserts++;
      stats.decoder.pagecache.count = pages_.count;
      pagelisthi->refcount--;
    }
    pagelisthi->refcount++;
    pagelisthi->add(bb, bb->mfnhi_loc);
    logging::println(logging::TRACE, "Add bb {} ({}, {} bytes) to high page list {}: loc {}:{}", (void*)bb,
                     (void*)(Waddr)bb->rip, bb->bytes, (void*)pagelisthi, (void*)bb->mfnhi_loc.chunk,
                     bb->mfnhi_loc.index);
    pagelisthi->refcount--;
  }

  pagelist->refcount--;

  logging::println(logging::DEBUG, "=====================================================================");
  logging::println(logging::DEBUG, "{}", *bb);
  logging::println(logging::DEBUG, "End of basic block: rip {} -> taken rip 0x{}, not taken rip 0x{}",
                   (void*)(Waddr)trans.bb.rip, (void*)(Waddr)trans.bb.rip_taken, (void*)(Waddr)trans.bb.rip_not_taken);

  bb->release();

  return bb;
}


#define MAX_RIP 0xffffffffffffffffULL

//
// Translate one basic block, just to get the uops: do not
// allocate an entry in the BB cache and do not add the BB
// to any lists. Instead, just copy the raw translated
// uops to targetbb and return.
//
// This function does not allocate any memory.
//
void BasicBlockCache::translate_in_place(BasicBlock& targetbb, Context& ctx, Waddr rip) {
  Options& config = ctx.machine_impl->config;

  if unlikely ((rip == config.log.start_log_at_rip) && (rip != MAX_RIP)) {
    config.log.start_log_at_iteration = 0;
    logenable = 1;
  }

  RIPVirtPhys rvp = rip;
  rvp.update(ctx);

  byte insnbuf[MAX_BB_BYTES];

  TraceDecoder trans(rvp);
  trans.fillbuf(ctx, insnbuf, sizeof(insnbuf));

  if (log_code_page_ops) {
    logging::println(logging::INFO, "Translating {} ({} bytes valid) at {} cycles, {} commits", W64(rvp),
                     trans.valid_byte_count, ctx.machine_impl->sim_cycle, ctx.machine_impl->total_user_insns_committed);
  }

  for (;;) {
    if (!trans.translate())
      break;
  }

  trans.bb.hitcount = 0;
  trans.bb.predcount = 0;

  memcpy(&targetbb, &trans.bb, sizeof(BasicBlockBase));
  memcpy(&targetbb.transops, &trans.bb.transops, trans.bb.count * sizeof(TransOp));

  logging::println(logging::DEBUG, "=====================================================================");
  logging::println(logging::DEBUG, "{}", targetbb);
  logging::println(logging::DEBUG, "End of basic block: rip {} -> taken rip 0x{}, not taken rip 0x{}",
                   (void*)(Waddr)trans.bb.rip, (void*)(Waddr)trans.bb.rip_taken, (void*)(Waddr)trans.bb.rip_not_taken);
}

BasicBlock* BasicBlockCache::translate_and_clone(Context& ctx, Waddr rip) {
  Options& config = ctx.machine_impl->config;

  if unlikely ((rip == config.log.start_log_at_rip) && (rip != MAX_RIP)) {
    config.log.start_log_at_iteration = 0;
    logenable = 1;
  }

  RIPVirtPhys rvp = rip;
  rvp.update(ctx);

  byte insnbuf[MAX_BB_BYTES];

  TraceDecoder trans(rvp);
  trans.fillbuf(ctx, insnbuf, sizeof(insnbuf));

  if (log_code_page_ops) {
    logging::println(logging::INFO, "Translating {} ({} bytes valid) at {} cycles, {} commits", W64(rvp),
                     trans.valid_byte_count, ctx.machine_impl->sim_cycle, ctx.machine_impl->total_user_insns_committed);
  }

  for (;;) {
    if (!trans.translate())
      break;
  }

  BasicBlock* bb = trans.bb.clone();

  return bb;
}

#undef MAX_RIP

} // namespace x86sim

auto std::formatter<x86sim::BasicBlockCache>::format(x86sim::BasicBlockCache& bbc, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  auto bblist = bbc.getentries();

  for (const auto* ptr : bblist) {
    const x86sim::BasicBlock& bb = *ptr;

    out = std::format_to(out,
                         "  {}: {:>4}t {:>3}ld {:>3}st {:>3}u {:>10}h {:>10}pr {:>10}uu : taken 0x{:012x}, seq {:012x}",
                         bb.rip, bb.tagcount, bb.memcount - bb.storecount, bb.storecount, bb.user_insn_count,
                         bb.hitcount, bb.predcount, bb.hitcount * bb.tagcount, bb.rip_taken, bb.rip_not_taken);
    //floatstring(100.0 * (double)bb.predcount / (double)bb.hitcount, 10, 2), "%pr "
    //floatstring(100.0 * percent_of_total_uops, 6, 2), "%uops",
    if (bb.rip_taken == bb.rip.rip)
      out = std::format_to(out, " [loop]");
    if (bb.repblock)
      out = std::format_to(out, " [repblock]");
    out = std::format_to(out, "\n");
  }

  return out;
}

namespace x86sim {

void shutdown_decode() {
  if (bbcache_dump_file) {
    std::fclose(bbcache_dump_file);
    bbcache_dump_file = nullptr;
  }
}

} // namespace x86sim
