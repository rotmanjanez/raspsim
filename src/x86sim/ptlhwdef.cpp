//
// PTLsim: Cycle Accurate x86-64 Simulator
// Hardware Definitions
//
// Copyright 1999-2008 Matt T. Yourst <yourst@yourst.com>
//

#include "ptlhwdef.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <format>
#include "globals.h"
#include <string>
#include "typedefs.h"

namespace x86sim {

//
// Micro-operation (uop) definitions
//

const OpcodeInfo opinfo[OP_MAX_OPCODE] = {
    // name, opclass, latency, fu
    {"nop", OPCLASS_LOGIC, opNOSIZE},
    {"mov", OPCLASS_LOGIC, opAB | ccB}, // move or merge
    // Logical
    {"and", OPCLASS_LOGIC, opAB},
    {"andnot", OPCLASS_LOGIC, opAB},
    {"xor", OPCLASS_LOGIC, opAB},
    {"or", OPCLASS_LOGIC, opAB},
    {"nand", OPCLASS_LOGIC, opAB},
    {"ornot", OPCLASS_LOGIC, opAB},
    {"eqv", OPCLASS_LOGIC, opAB},
    {"nor", OPCLASS_LOGIC, opAB},
    // Mask, insert or extract bytes
    {"maskb", OPCLASS_SIMPLE_SHIFT,
     opABC}, // mask rd = ra,rb,[ds,ms,mc], bytes only; rcimm (8 bits, but really 18 bits)
    // Add and subtract
    {"add", OPCLASS_ADDSUB, opABC | ccC}, // ra + rb
    {"sub", OPCLASS_ADDSUB, opABC | ccC}, // ra - rb
    {"adda", OPCLASS_ADDSHIFT, opABC},    // ra + rb + rc
    {"suba", OPCLASS_ADDSHIFT, opABC},    // ra - rb + rc
    {"addm", OPCLASS_ADDSUB, opABC},      // lowbits(ra + rb, m)
    {"subm", OPCLASS_ADDSUB, opABC},      // lowbits(ra - rb, m)
    // Condition code logical ops
    {"andcc", OPCLASS_FLAGS, opAB | ccAB | opNOSIZE},
    {"orcc", OPCLASS_FLAGS, opAB | ccAB | opNOSIZE},
    {"xorcc", OPCLASS_FLAGS, opAB | ccAB | opNOSIZE},
    {"ornotcc", OPCLASS_FLAGS, opAB | ccAB | opNOSIZE},
    // Condition code movement and merging
    {"movccr", OPCLASS_FLAGS, opB | ccB | opNOSIZE},
    {"movrcc", OPCLASS_FLAGS, opB | opNOSIZE},
    {"collcc", OPCLASS_FLAGS, opABC | ccABC | opNOSIZE},
    // Simple shifting (restricted to small immediate 1..8)
    {"shls", OPCLASS_SIMPLE_SHIFT, opAB}, // rb imm limited to 0-8
    {"shrs", OPCLASS_SIMPLE_SHIFT, opAB}, // rb imm limited to 0-8
    {"bswap", OPCLASS_LOGIC, opAB},       // byte swap rb
    {"sars", OPCLASS_SIMPLE_SHIFT, opAB}, // rb imm limited to 0-8
    // Bit testing
    {"bt", OPCLASS_LOGIC, opAB},
    {"bts", OPCLASS_LOGIC, opAB},
    {"btr", OPCLASS_LOGIC, opAB},
    {"btc", OPCLASS_LOGIC, opAB},
    // Set and select
    {"set", OPCLASS_SELECT, opABC | ccAB},     // rd = rc <- (eval(ra,rb) ? 1 : 0)
    {"set.sub", OPCLASS_SELECT, opABC},        // rd = rc <- (eval(ra-rb) ? 1 : 0)
    {"set.and", OPCLASS_SELECT, opABC},        // rd = rc <- (eval(ra&rb) ? 1 : 0)
    {"sel", OPCLASS_SELECT, opABC | ccABC},    // rd = falsereg,truereg,condreg
    {"sel.cmp", OPCLASS_SELECT, opABC | ccAB}, // rd = falsereg,truereg,intreg
    // Branches
    {"br", OPCLASS_COND_BRANCH, opAB | ccAB | opNOSIZE}, // branch (rcimm: 32 to 53-bit target info)
    {"br.sub", OPCLASS_COND_BRANCH,
     opAB}, // compare and branch ("cmp" form: subtract) (rcimm: 32 to 53-bit target info)
    {"br.and", OPCLASS_COND_BRANCH, opAB},    // compare and branch ("test" form: and) (rcimm: 32 to 53-bit target info)
    {"jmp", OPCLASS_INDIR_BRANCH, opA},       // indirect user branch (rcimm: 32 to 53-bit target info)
    {"bru", OPCLASS_UNCOND_BRANCH, opNOSIZE}, // unconditional branch (rcimm: 32 to 53-bit target info)
    {"jmpp", OPCLASS_INDIR_BRANCH | OPCLASS_BARRIER,
     opA}, // indirect branch within PTL (rcimm: 32 to 53-bit target info)
    {"brp", OPCLASS_UNCOND_BRANCH | OPCLASS_BARRIER,
     opNOSIZE}, // unconditional branch (PTL only) (rcimm: 32 to 53-bit target info)
    // Checks
    {"chk", OPCLASS_CHECK,
     opABC | ccAB | opNOSIZE}, // check condition and rollback if false (uses cond codes); (rcimm: 8-bit exception type)
    {"chk.sub", OPCLASS_CHECK, opABC}, // check ("cmp" form: subtract)
    {"chk.and", OPCLASS_CHECK, opABC}, // check ("test" form: and)
    // Loads and stores
    {"ld", OPCLASS_LOAD, opABC},        // load zero extended
    {"ldx", OPCLASS_LOAD, opABC},       // load sign extended
    {"ld.pre", OPCLASS_PREFETCH, opAB}, // prefetch
    {"ld.a16", OPCLASS_LOAD, opABC},    // load zero extended, 16-byte aligned
    {"st", OPCLASS_STORE, opABC},       // store
    {"st.a16", OPCLASS_STORE, opABC},   // store, 16-byte aligned
    {"mf", OPCLASS_FENCE, opNOSIZE},    // memory fence (extshift holds type: 01 = st, 10 = ld, 11 = ld.st)
    // Shifts, rotates and complex masking
    {"shl", OPCLASS_SHIFTROT, opABC | ccC},
    {"shr", OPCLASS_SHIFTROT, opABC | ccC},
    {"mask", OPCLASS_SHIFTROT, opAB}, // mask rd = ra,rb,[ds,ms,mc]: (rcimm: 18 bits)
    {"sar", OPCLASS_SHIFTROT, opABC | ccC},
    {"rotl", OPCLASS_SHIFTROT, opABC | ccC},
    {"rotr", OPCLASS_SHIFTROT, opABC | ccC},
    {"rotcl", OPCLASS_SHIFTROT, opABC | ccC},
    {"rotcr", OPCLASS_SHIFTROT, opABC | ccC},
    // Multiplication
    {"mull", OPCLASS_MULTIPLY, opAB},
    {"mulh", OPCLASS_MULTIPLY, opAB},
    {"mulhu", OPCLASS_MULTIPLY, opAB},
    {"mulhl", OPCLASS_MULTIPLY, opAB},
    // Bit scans
    {"ctz", OPCLASS_BITSCAN, opB},
    {"clz", OPCLASS_BITSCAN, opB},
    {"ctpop", OPCLASS_BITSCAN, opB},
    {"permb", OPCLASS_SHIFTROT, opABC},
    // Integer divide and remainder step
    {"div", OPCLASS_MULTIPLY, opABC},  // unsigned divide
    {"rem", OPCLASS_MULTIPLY, opABC},  // unsigned divide
    {"divs", OPCLASS_MULTIPLY, opABC}, // signed divide
    {"rems", OPCLASS_MULTIPLY, opABC}, // signed divide
    // Minimum and maximum
    {"min", OPCLASS_ADDSUB, opAB},   // min(ra, rb)
    {"max", OPCLASS_ADDSUB, opAB},   // max(ra, rb)
    {"min.s", OPCLASS_ADDSUB, opAB}, // min(ra, rb) (ra and rb are signed types)
    {"max.s", OPCLASS_ADDSUB, opAB}, // max(ra, rb) (ra and rb are signed types)
    // Floating point
    // uop.size bits have following meaning:
    // 00 = single precision, scalar (preserve high 32 bits of ra)
    // 01 = single precision, packed (two 32-bit floats)
    // 1x = double precision, scalar or packed (use two uops to process 128-bit xmm)
    {"fadd", OPCLASS_FP_ALU, opAB},
    {"fsub", OPCLASS_FP_ALU, opAB},
    {"fmul", OPCLASS_FP_ALU, opAB},
    {"fmadd", OPCLASS_FP_ALU, opABC},
    {"fmsub", OPCLASS_FP_ALU, opABC},
    {"fmsubr", OPCLASS_FP_ALU, opABC},
    {"fdiv", OPCLASS_FP_DIVSQRT, opAB},
    {"fsqrt", OPCLASS_FP_DIVSQRT, opAB},
    {"frcp", OPCLASS_FP_DIVSQRT, opAB},
    {"fsqrt", OPCLASS_FP_DIVSQRT, opAB},
    {"fmin", OPCLASS_FP_COMPARE, opAB},
    {"fmax", OPCLASS_FP_COMPARE, opAB},
    {"fcmp", OPCLASS_FP_COMPARE, opAB},
    // For fcmpcc, uop.size bits have following meaning:
    // 00 = single precision ordered compare
    // 01 = single precision unordered compare
    // 10 = double precision ordered compare
    // 11 = double precision unordered compare
    {"fcmpcc", OPCLASS_FP_COMPARE, opAB},
    // and/andn/or/xor are done using integer uops
    // For these conversions, uop.size bits select truncation mode:
    // x0 = normal IEEE-style rounding
    // x1 = truncate to zero
    {"fcvt.i2s.ins", OPCLASS_FP_CONVERTI2F,
     opAB}, // one W32s <rb> to single, insert into low 32 bits of <ra> (for cvtsi2ss)
    {"fcvt.i2s.p", OPCLASS_FP_CONVERTI2F, opB}, // pair of W32s <rb> to pair of singles <rd> (for cvtdq2ps, cvtpi2ps)
    {"fcvt.i2d.lo", OPCLASS_FP_CONVERTI2F,
     opB}, // low W32s in <rb> to double in <rd> (for cvtdq2pd part 1, cvtpi2pd part 1, cvtsi2sd)
    {"fcvt.i2d.hi", OPCLASS_FP_CONVERTI2F,
     opB}, // high W32s in <rb> to double in <rd> (for cvtdq2pd part 2, cvtpi2pd part 2)
    {"fcvt.q2s.ins", OPCLASS_FP_CONVERTI2F,
     opAB}, // one W64s <rb> to single, insert into low 32 bits of <ra> (for cvtsi2ss with REX.mode64 prefix)
    {"fcvt.q2d", OPCLASS_FP_CONVERTI2F,
     opAB}, // one W64s <rb> to double in <rd>, ignore <ra> (for cvtsi2sd with REX.mode64 prefix)
    {"fcvt.s2i", OPCLASS_FP_CONVERTF2I, opB}, // one single <rb> to W32s in <rd> (for cvtss2si, cvttss2si)
    {"fcvt.s2q", OPCLASS_FP_CONVERTF2I,
     opB}, // one single <rb> to W64s in <rd> (for cvtss2si, cvttss2si with REX.mode64 prefix)
    {"fcvt.s2i.p", OPCLASS_FP_CONVERTF2I,
     opB}, // pair of singles in <rb> to pair of W32s in <rd> (for cvtps2pi, cvttps2pi, cvtps2dq, cvttps2dq)
    {"fcvt.d2i", OPCLASS_FP_CONVERTF2I, opB}, // one double <rb> to W32s in <rd> (for cvtsd2si, cvttsd2si)
    {"fcvt.d2q", OPCLASS_FP_CONVERTF2I, opB}, // one double <rb> to W64s in <rd> (for cvtsd2si with REX.mode64 prefix)
    {"fcvt.d2i.p", OPCLASS_FP_CONVERTF2I,
     opAB}, // pair of doubles in <ra> (high), <rb> (low) to pair of W32s in <rd> (for cvtpd2pi, cvttpd2pi, cvtpd2dq, cvttpd2dq), clear high 64 bits of dest xmm
    {"fcvt.d2s.ins", OPCLASS_FP_CONVERTFP,
     opAB}, // double in <rb> to single, insert into low 32 bits of <ra> (for cvtsd2ss)
    {"fcvt.d2s.p", OPCLASS_FP_CONVERTFP,
     opAB}, // pair of doubles in <ra> (high), <rb> (low) to pair of singles in <rd> (for cvtpd2ps)
    {"fcvt.s2d.lo", OPCLASS_FP_CONVERTFP, opB}, // low single in <rb> to double in <rd> (for cvtps2pd, part 1, cvtss2sd)
    {"fcvt.s2d.hi", OPCLASS_FP_CONVERTFP, opB}, // high single in <rb> to double in <rd> (for cvtps2pd, part 2)
    // Vector integer uops
    // uop.size defines element size: 00 = byte, 01 = W16, 10 = W32, 11 = W64 (i.e. same as normal ALU uops)
    {"vadd", OPCLASS_VEC_ALU, opAB},    // vector add with wraparound
    {"vsub", OPCLASS_VEC_ALU, opAB},    // vector sub with wraparound
    {"vadd.us", OPCLASS_VEC_ALU, opAB}, // vector add with unsigned saturation
    {"vsub.us", OPCLASS_VEC_ALU, opAB}, // vector sub with unsigned saturation
    {"vadd.ss", OPCLASS_VEC_ALU, opAB}, // vector add with signed saturation
    {"vsub.ss", OPCLASS_VEC_ALU, opAB}, // vector sub with signed saturation
    {"vshl", OPCLASS_VEC_ALU, opAB},    // vector shift left
    {"vshr", OPCLASS_VEC_ALU, opAB},    // vector shift right
    {"vbt", OPCLASS_VEC_ALU, opAB}, // vector bit test (pack bit <rb> of each element in <ra> into low N bits of output)
    {"vsar", OPCLASS_VEC_ALU, opAB}, // vector shift right arithmetic (sign extend)
    {"vavg", OPCLASS_VEC_ALU, opAB}, // vector average ((<ra> + <rb> + 1) >> 1)
    {"vcmp", OPCLASS_VEC_ALU,
     opAB}, // vector compare (uop.cond specifies compare type; result is all 1's for true, or all 0's for false in each element)
    {"vmin", OPCLASS_VEC_ALU, opAB},     // vector minimum
    {"vmax", OPCLASS_VEC_ALU, opAB},     // vector maximum
    {"vmin.s", OPCLASS_VEC_ALU, opAB},   // vector signed minimum
    {"vmax.s", OPCLASS_VEC_ALU, opAB},   // vector signed maximum
    {"vmull", OPCLASS_VEC_ALU, opAB},    // multiply and keep low bits
    {"vmulh", OPCLASS_VEC_ALU, opAB},    // multiply and keep high bits
    {"vmulhu", OPCLASS_VEC_ALU, opAB},   // multiply and keep high bits (unsigned)
    {"vmaddp", OPCLASS_VEC_ALU, opAB},   // multiply and add adjacent pairs (signed)
    {"vsad", OPCLASS_VEC_ALU, opAB},     // sum of absolute differences
    {"vpack.us", OPCLASS_VEC_ALU, opAB}, // pack larger to smaller (unsigned saturation)
    {"vpack.ss", OPCLASS_VEC_ALU, opAB}, // pack larger to smaller (signed saturation)
};

const char* exception_names[EXCEPTION_COUNT] = {
    // 0123456789abcdef
    "NoException", "Propagate",     "BranchMiss",     "Unaligned",   "PageRead",  "PageWrite",
    "PageExec",    "StStAlias",     "LdStAlias",      "CheckFailed", "SkipBlock", "LFRQFull",
    "Float",       "FloatNotAvail", "DivideOverflow", "InvalidAddr",
};

const char* x86_exception_names[256] = {
    "divide",         "debug",         "nmi",          "breakpoint",     "overflow",    "bounds",
    "invalid opcode", "fpu not avail", "double fault", "coproc overrun", "invalid tss", "seg not present",
    "stack fault",    "gp fault",      "page fault",   "spurious int",   "fpu",         "unaligned",
    "machine check",  "sse",           "int14h",       "int15h",         "int16h",      "int17h",
    "int18h",         "int19h",        "int1Ah",       "int1Bh",         "int1Ch",      "int1Dh",
    "int1Eh",         "int1Fh",        "int20h",       "int21h",         "int22h",      "int23h",
    "int24h",         "int25h",        "int26h",       "int27h",         "int28h",      "int29h",
    "int2Ah",         "int2Bh",        "int2Ch",       "int2Dh",         "int2Eh",      "int2Fh",
    "int30h",         "int31h",        "int32h",       "int33h",         "int34h",      "int35h",
    "int36h",         "int37h",        "int38h",       "int39h",         "int3Ah",      "int3Bh",
    "int3Ch",         "int3Dh",        "int3Eh",       "int3Fh",         "int40h",      "int41h",
    "int42h",         "int43h",        "int44h",       "int45h",         "int46h",      "int47h",
    "int48h",         "int49h",        "int4Ah",       "int4Bh",         "int4Ch",      "int4Dh",
    "int4Eh",         "int4Fh",        "int50h",       "int51h",         "int52h",      "int53h",
    "int54h",         "int55h",        "int56h",       "int57h",         "int58h",      "int59h",
    "int5Ah",         "int5Bh",        "int5Ch",       "int5Dh",         "int5Eh",      "int5Fh",
    "int60h",         "int61h",        "int62h",       "int63h",         "int64h",      "int65h",
    "int66h",         "int67h",        "int68h",       "int69h",         "int6Ah",      "int6Bh",
    "int6Ch",         "int6Dh",        "int6Eh",       "int6Fh",         "int70h",      "int71h",
    "int72h",         "int73h",        "int74h",       "int75h",         "int76h",      "int77h",
    "int78h",         "int79h",        "int7Ah",       "int7Bh",         "int7Ch",      "int7Dh",
    "int7Eh",         "int7Fh",        "int80h",       "int81h",         "int82h",      "int83h",
    "int84h",         "int85h",        "int86h",       "int87h",         "int88h",      "int89h",
    "int8Ah",         "int8Bh",        "int8Ch",       "int8Dh",         "int8Eh",      "int8Fh",
    "int90h",         "int91h",        "int92h",       "int93h",         "int94h",      "int95h",
    "int96h",         "int97h",        "int98h",       "int99h",         "int9Ah",      "int9Bh",
    "int9Ch",         "int9Dh",        "int9Eh",       "int9Fh",         "intA0h",      "intA1h",
    "intA2h",         "intA3h",        "intA4h",       "intA5h",         "intA6h",      "intA7h",
    "intA8h",         "intA9h",        "intAAh",       "intABh",         "intACh",      "intADh",
    "intAEh",         "intAFh",        "intB0h",       "intB1h",         "intB2h",      "intB3h",
    "intB4h",         "intB5h",        "intB6h",       "intB7h",         "intB8h",      "intB9h",
    "intBAh",         "intBBh",        "intBCh",       "intBDh",         "intBEh",      "intBFh",
    "intC0h",         "intC1h",        "intC2h",       "intC3h",         "intC4h",      "intC5h",
    "intC6h",         "intC7h",        "intC8h",       "intC9h",         "intCAh",      "intCBh",
    "intCCh",         "intCDh",        "intCEh",       "intCFh",         "intD0h",      "intD1h",
    "intD2h",         "intD3h",        "intD4h",       "intD5h",         "intD6h",      "intD7h",
    "intD8h",         "intD9h",        "intDAh",       "intDBh",         "intDCh",      "intDDh",
    "intDEh",         "intDFh",        "intE0h",       "intE1h",         "intE2h",      "intE3h",
    "intE4h",         "intE5h",        "intE6h",       "intE7h",         "intE8h",      "intE9h",
    "intEAh",         "intEBh",        "intECh",       "intEDh",         "intEEh",      "intEFh",
    "intF0h",         "intF1h",        "intF2h",       "intF3h",         "intF4h",      "intF5h",
    "intF6h",         "intF7h",        "intF8h",       "intF9h",         "intFAh",      "intFBh",
    "intFCh",         "intFDh",        "intFEh",       "intFFh"};

const char* arch_reg_names[TRANSREG_COUNT] = {
    // Integer registers
    "rax",
    "rcx",
    "rdx",
    "rbx",
    "rsp",
    "rbp",
    "rsi",
    "rdi",
    "r8",
    "r9",
    "r10",
    "r11",
    "r12",
    "r13",
    "r14",
    "r15",
    // SSE registers
    "xmml0",
    "xmmh0",
    "xmml1",
    "xmmh1",
    "xmml2",
    "xmmh2",
    "xmml3",
    "xmmh3",
    "xmml4",
    "xmmh4",
    "xmml5",
    "xmmh5",
    "xmml6",
    "xmmh6",
    "xmml7",
    "xmmh7",
    "xmml8",
    "xmmh8",
    "xmml9",
    "xmmh9",
    "xmml10",
    "xmmh10",
    "xmml11",
    "xmmh11",
    "xmml12",
    "xmmh12",
    "xmml13",
    "xmmh13",
    "xmml14",
    "xmmh14",
    "xmml15",
    "xmmh15",
    // x87 FP/MMX
    "fptos",
    "fpsw",
    "fptags",
    "fpstack",
    "msr",
    "dlptr",
    "trace",
    "ctx",
    // Special
    "rip",
    "flags",
    "dlend",
    "selfrip",
    "nextrip",
    "ar1",
    "ar2",
    "zero",
    // The following are ONLY used during the translation and renaming process:
    "tr0",
    "tr1",
    "tr2",
    "tr3",
    "tr4",
    "tr5",
    "tr6",
    "tr7",
    "zf",
    "cf",
    "of",
    "imm",
    "mem",
    "tr8",
    "tr9",
    "tr10",
};

void Context::fxsave(FXSAVEStruct& state) {
  state.cw = static_cast<W16>(fpcw);
  // clear everything but 4 FP status flag bits (c3/c2/c1/c0):
  state.sw = commitarf[REG_fpsw] & ((0x7 << 8) | (1 << 14));
  int tos = commitarf[REG_fptos] >> 3;
  assert(inrange(tos, 0, 7));
  state.sw.tos = tos;
  state.tw = 0;

  // Prepare tag word (special format for FXSAVE)
  foreach (i, 8)
    state.tw |= (bit(commitarf[REG_fptags], i * 8) << i);

  // Prepare actual registers
  foreach (i, 8)
    x87_fp_64bit_to_80bit(&state.fpregs[i].reg, fpstack[lowbits(tos + i, 3)]);

  state.fop = 0;

  if (use64) {
    state.use64.rip = 0;
    state.use64.rdp = 0;
  } else {
    state.use32.eip = 0;
    state.use32.cs = 0;
    state.use32.dp = 0;
    state.use32.ds = 0;
  }

  state.mxcsr = mxcsr;
  state.mxcsr_mask = 0x0000ffff; // all MXCSR features supported

  foreach (i, (use64) ? 16 : 8) {
    state.xmmregs[i].lo = commitarf[REG_xmml0 + i * 2];
    state.xmmregs[i].hi = commitarf[REG_xmmh0 + i * 2];
  }
}

void Context::fxrstor(const FXSAVEStruct& state) {
  commitarf[REG_fptos] = state.sw.tos * 8;
  commitarf[REG_fpsw] = state.sw;
  fpcw = state.cw;

  commitarf[REG_fptags] = 0;
  foreach (i, 8) {
    // FXSAVE struct uses an abbreviated tag word with 8 bits (0 = empty, 1 = used)
    int used = bit(state.tw, i);
    commitarf[REG_fptags] |= ((W64)used) << i * 8;
  }

  // x86 FSAVE state is in order of stack rather than physical registers:
  foreach (i, 8) {
    fpstack[lowbits(state.sw.tos + i, 3)] =
        x87_fp_80bit_to_64bit(&state.fpregs[i].reg, X87ControlWord(static_cast<W16>(fpcw)).rc);
  }

  mxcsr = state.mxcsr & state.mxcsr_mask;

  foreach (i, (use64) ? 16 : 8) {
    commitarf[REG_xmml0 + i * 2] = state.xmmregs[i].lo;
    commitarf[REG_xmmh0 + i * 2] = state.xmmregs[i].hi;
  }
}

/*
 * Convert a condition code (as in jump, setcc, cmovcc, etc) to
 * the one or two architectural registers last updated with the
 * flags that uop will test.
 */
const CondCodeToFlagRegs cond_code_to_flag_regs[16] = {
    {0, REG_of, REG_of}, // of:               jo          (rb only)
    {0, REG_of, REG_of}, // !of:              jno         (rb only)
    {0, REG_cf, REG_cf}, // cf:               jb jc jnae  (rb only)
    {0, REG_cf, REG_cf}, // !cf:              jnb jnc jae (rb only)
    {0, REG_zf, REG_zf}, // zf:               jz je       (ra only)
    {0, REG_zf, REG_zf}, // !zf:              jnz jne     (ra only)
    {1, REG_zf, REG_cf}, // cf|zf:            jbe jna
    {1, REG_zf, REG_cf}, // !cf & !zf:        jnbe ja
    {0, REG_zf, REG_zf}, // sf:               js          (ra only)
    {0, REG_zf, REG_zf}, // !sf:              jns         (ra only)
    {0, REG_zf, REG_zf}, // pf:               jp jpe      (ra only)
    {0, REG_zf, REG_zf}, // !pf:              jnp jpo     (ra only)
    {1, REG_zf, REG_of}, // sf != of:         jl jnge (*)
    {1, REG_zf, REG_of}, // sf == of:         jnl jge (*)
    {1, REG_zf, REG_of}, // zf | (sf != of):  jle jng (*)
    {1, REG_zf, REG_of}, // !zf & (sf == of): jnle jg (*)
                         //
                         // (*) Technically three flags are involved in the comparison here,
                         // however as pursuant to the ZAPS trick, zf/af/pf/sf are always
                         // either all written together or not written at all. Hence the
                         // last writer of SF will also deliver ZF in the same result.
                         //
};

const char* cond_code_names[16] = {"o", "no", "c", "nc", "e", "ne", "be", "nbe",
                                   "s", "ns", "p", "np", "l", "nl", "le", "nle"};
const char* x86_flag_names[32] = {"c",  "X",     "p",     "W",  "a",  "U",  "z",  "s",  "t",   "i",   "d",
                                  "o",  "iopl0", "iopl1", "nt", "B",  "rf", "vm", "ac", "vif", "vip", "id",
                                  "22", "23",    "24",    "25", "26", "27", "28", "29", "30",  "31"};

const char* setflag_names[SETFLAG_COUNT] = {"z", "c", "o"};
const W16 setflags_to_x86_flags[1 << 3] = {
    0 | 0 | 0,                     // 000 = n/a
    0 | 0 | FLAG_ZAPS,             // 001 = Z
    0 | FLAG_CF | 0,               // 010 =  C
    0 | FLAG_CF | FLAG_ZAPS,       // 011 = ZC
    FLAG_OF | 0 | 0,               // 100 =   O
    FLAG_OF | 0 | FLAG_ZAPS,       // 101 = Z O
    FLAG_OF | FLAG_CF | 0,         // 110 =  CO
    FLAG_OF | FLAG_CF | FLAG_ZAPS, // 111 = ZCO
};

void BasicBlock::reset() {
  setzero(*((BasicBlockBase*)this));
  hashlink.reset();
  mfnlo_loc.reset();
  mfnhi_loc.reset();
  type = BB_TYPE_COND;
}

void BasicBlock::reset(const RIPVirtPhys& rip) {
  reset();
  this->rip = rip;
  rip_taken = rip;
  rip_not_taken = rip;
}

//
// This is explicitly defined instead of just using a
// destructor since we do some fancy dynamic resizing
// in the clone() method that c++ will croak on.
//
// Once you call this, the basic block is *gone* and
// cannot be accessed ever again, even if it is still
// in scope. Don't call this with non-cloned() blocks.
//
void BasicBlock::free() {
  if (synthops)
    delete[] synthops;
  synthops = null;
  ::free(this);
}

BasicBlock* BasicBlock::clone() {
  BasicBlock* bb = (BasicBlock*)malloc(sizeof(BasicBlockBase) + (count * sizeof(TransOp)));

  memcpy(bb, this, sizeof(BasicBlockBase));

  bb->synthops = null;
  // hashlink, mfnlo_loc, mfnhi_loc are always updated after cloning
  bb->hashlink.reset();
  bb->use(0);

  foreach (i, count)
    bb->transops[i] = this->transops[i];
  return bb;
}

const char* bb_type_names[BB_TYPE_COUNT] = {"br", "bru", "jmp", "brp"};

char* regname(int r) {
  static char temp[16];
  assert(r >= 0);
  assert(r < 256);
  std::snprintf(temp, sizeof(temp), "r%d", r);
  return temp;
}

std::string nameof(const TransOpBase& uop) {
  static const char* size_names[4] = {"b", "w", "d", ""};
  static const char* fptype_names[4] = {"ss", "ps", "sd", "pd"};
  static const char* mask_exttype[4] = {"", "zxt", "sxt", "???"};

  int op = uop.opcode;

  bool ld = isload(op);
  bool st = isstore(op);
  bool fp = (isclass(op, OPCLASS_FP_ALU));

  std::string result = nameof(op);

  if ((op != OP_maskb) & (op != OP_mask))
    result += (fp ? fptype_names[uop.size] : size_names[uop.size]);
  else
    result += std::format(".{}", mask_exttype[uop.cond]);

  if (isclass(op, OPCLASS_USECOND))
    result += std::format(".{}", cond_code_names[uop.cond]);

  if (ld | st) {
    result += ((uop.cond == LDST_ALIGN_LO) ? ".low" : (uop.cond == LDST_ALIGN_HI) ? ".high" : "");
    if (uop.cachelevel > 0)
      result += std::format(".L{}", (char)('1' + uop.cachelevel));
  }

  if (uop.internal)
    result += ".p";

  return result;
}

} // namespace x86sim

auto std::formatter<x86sim::UserContext>::format(const x86sim::UserContext& arf, std::format_context& ctx) const {
  using namespace x86sim;
  static const int width = 4;
  auto out = ctx.out();
  foreach (i, ARCHREG_COUNT) {
    out = std::format_to(out, "  {:<6} 0x{:016x}  ", arch_reg_names[i], arf[i]);
    if ((i % width) == (width - 1))
      out = std::format_to(out, "\n");
  }
  return out;
}

namespace x86sim {

std::string format_value_and_flags(W64 value, W16 flags) {
  std::string flagstr;
  if (flags & FLAG_ZF)
    flagstr += 'z';
  if (flags & FLAG_PF)
    flagstr += 'p';
  if (flags & FLAG_SF)
    flagstr += 's';
  if (flags & FLAG_CF)
    flagstr += 'c';
  if (flags & FLAG_OF)
    flagstr += 'o';

  std::string result;
  if (flags & FLAG_INV) {
    W32 exc_idx = LO32(value);
    const char* exc_name = (exc_idx < 256) ? x86_exception_names[exc_idx] : "unknown";
    result = std::format(" < {:<14} >", exc_name);
  } else {
    result = std::format(" 0x{:016x}", value);
  }
  result += std::format("|{:<5}", flagstr);
  return result;
}

} // namespace x86sim
