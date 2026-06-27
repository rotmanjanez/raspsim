//
// PTLsim: Cycle Accurate x86-64 Simulator
// Interface to uop implementations
//
// Copyright 2000-2008 Matt T. Yourst <yourst@yourst.com>
//

#include "globals.h"
#include "ptlhwdef.h"
#include "ptlsim.h"
#include "typedefs.h"
#include "x86sim/logging.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>


namespace x86sim {

//
// Flags generation (all but CF and OF)
//
template<typename T>
inline byte x86_genflags(T r) {
  using U = std::make_unsigned_t<T>;

  const U value = static_cast<U>(r);
  byte flags = 0;

  if ((value >> ((sizeof(T) * 8) - 1)) & 1)
    flags |= FLAG_SF;
  if (value == 0)
    flags |= FLAG_ZF;
  if ((std::popcount(static_cast<unsigned>(value & 0xff)) & 1) == 0)
    flags |= FLAG_PF;

  return flags;
}

template byte x86_genflags<byte>(byte r);
template byte x86_genflags<W16>(W16 r);
template byte x86_genflags<W32>(W32 r);
template byte x86_genflags<W64>(W64 r);

//
// Flags format: OF - - - SF ZF - AF - PF - CF
//               11       7  6    4    2    0
//               rb       ra ra   ra   ra   rb
//

template<typename T>
inline W64 x86_merge(W64 rd, W64 ra) {
  union {
    W8 w8;
    W16 w16;
    W32 w32;
    W64 w64;
  } sizes;

  switch (sizeof(T)) {
  case 1:
    sizes.w64 = rd;
    sizes.w8 = ra;
    return sizes.w64;
  case 2:
    sizes.w64 = rd;
    sizes.w16 = ra;
    return sizes.w64;
  case 4:
    return LO32(ra);
  case 8:
    return ra;
  }

  return rd;
}

#define ZAPS SETFLAG_ZF
#define CF SETFLAG_CF
#define OF SETFLAG_OF

#define make_exp_aluop(name, expr)                                                                                     \
  template<typename T, int genflags>                                                                                   \
  struct name {                                                                                                        \
    T operator()(T ra, T rb, T rc, W16 raflags, W16 rbflags, W16 rcflags, byte& cf, byte& of) {                        \
      cf = 0;                                                                                                          \
      of = 0;                                                                                                          \
      W64 rd;                                                                                                          \
      expr;                                                                                                            \
      return rd;                                                                                                       \
    }                                                                                                                  \
  }

#define x86_aluop2_opcode_add std::plus
#define x86_aluop2_opcode_adc std::plus
#define x86_aluop2_opcode_sub std::minus
#define x86_aluop2_opcode_sbb std::minus

#define x86_aluop2_carry_add false
#define x86_aluop2_carry_adc true
#define x86_aluop2_carry_sub false
#define x86_aluop2_carry_sbb true

template<template<typename> class Operation, bool carry_in, bool input_flags, typename T, int genflags>
inline T x86_aluop2(T ra, T rb, W16 rcflags, byte& cf, byte& of) {
  using U = std::make_unsigned_t<T>;

  const U a = static_cast<U>(ra);
  const U b = static_cast<U>(rb);
  const U carry = (input_flags && carry_in && ((rcflags & FLAG_CF) != 0)) ? U(1) : U(0);
  constexpr bool subtract = std::is_same_v<Operation<U>, std::minus<U>>;
  constexpr U sign_bit = U(1) << (std::numeric_limits<U>::digits - 1);
  const U intermediate = U(Operation<U>{}(a, b));
  const U result = U(Operation<U>{}(intermediate, carry));
  const bool carry_out =
      subtract ? ((a < b) || (intermediate < carry)) : ((intermediate < a) || (result < intermediate));
  const U rhs_sign = subtract ? U(~b) : b;
  const bool overflow = (~(a ^ rhs_sign) & (a ^ result) & sign_bit) != 0;

  cf = (genflags & SETFLAG_CF) ? byte(carry_out) : 0;
  of = (genflags & SETFLAG_OF) ? byte(overflow) : 0;
  return static_cast<T>(result);
}

#define make_x86_aluop2(name, opcode, pretext)                                                                         \
  template<typename T, int genflags>                                                                                   \
  struct name {                                                                                                        \
    T operator()(T ra, T rb, T rc, W16 raflags, W16 rbflags, W16 rcflags, byte& cf, byte& of) {                        \
      return x86_aluop2<x86_aluop2_opcode_##opcode, x86_aluop2_carry_##opcode, (sizeof(pretext) > 1), T, genflags>(    \
          ra, rb, rcflags, cf, of);                                                                                    \
    }                                                                                                                  \
  }

UopResult uop_impl_nop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = 0;
  return out;
}

//
// 2-operand ALU operation
//
template<int ptlopcode, template<typename, int> class func, typename T, int genflags>
inline UopResult aluop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  byte cf = 0, of = 0;
  func<T, genflags> f;
  T rt = f(ra, rb, rc, raflags, rbflags, rcflags, cf, of);
  out.rddata = x86_merge<T>(ra, rt);
  out.rdflags = (of << 11) | cf | ((genflags & SETFLAG_ZF) ? x86_genflags<T>(rt) : 0);
  return out;
}

#define make_anyop_all_sizes(ptlopcode, mapname, opclass, nativeop, flagset)                                           \
  UopImpl mapname[4][2] = {{&opclass<ptlopcode, nativeop, W8, 0>, &opclass<ptlopcode, nativeop, W8, (flagset)>},       \
                           {&opclass<ptlopcode, nativeop, W16, 0>, &opclass<ptlopcode, nativeop, W16, (flagset)>},     \
                           {&opclass<ptlopcode, nativeop, W32, 0>, &opclass<ptlopcode, nativeop, W32, (flagset)>},     \
                           {&opclass<ptlopcode, nativeop, W64, 0>, &opclass<ptlopcode, nativeop, W64, (flagset)>}}

#define make_aluop_all_sizes(ptlopcode, mapname, nativeop, flagset)                                                    \
  make_anyop_all_sizes(ptlopcode, mapname, aluop, nativeop, flagset);

#define make_exp_aluop_all_sizes(name, exp, setflags)                                                                  \
  make_exp_aluop(exp_op_##name, (exp));                                                                                \
  make_aluop_all_sizes(OP_##name, implmap_##name, exp_op_##name, (setflags));

#define make_x86_aluop_all_sizes(name, opcode, setflags, pretext)                                                      \
  make_x86_aluop2(x86_op_##name, opcode, pretext);                                                                     \
  make_aluop_all_sizes(OP_##name, implmap_##name, x86_op_##name, (setflags));

#define PRETEXT_NO_FLAGS_IN ""
#define PRETEXT_ALL_FLAGS_IN "pushw %[rcflags]; popfw; "

template<typename T>
inline UopResult exp_op_mov(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = x86_merge<T>(ra, rb);
  out.rdflags = rbflags;
  return out;
}

UopImpl implmap_mov[4] = {&exp_op_mov<W8>, &exp_op_mov<W16>, &exp_op_mov<W32>, &exp_op_mov<W64>};

// make_exp_aluop_all_sizes(mov, (rd = (rb)), 0);
make_exp_aluop_all_sizes(and, (rd = (ra & rb)), ZAPS);
make_exp_aluop_all_sizes(or, (rd = (ra | rb)), ZAPS);
make_exp_aluop_all_sizes(xor, (rd = (ra ^ rb)), ZAPS);
make_exp_aluop_all_sizes(andnot, (rd = ((~ra) & rb)), ZAPS);
make_exp_aluop_all_sizes(ornot, (rd = ((~ra) | rb)), ZAPS);
make_exp_aluop_all_sizes(nand, (rd = (~(ra & rb))), ZAPS);
make_exp_aluop_all_sizes(nor, (rd = (~(ra | rb))), ZAPS);
make_exp_aluop_all_sizes(eqv, (rd = (~(ra ^ rb))), ZAPS);
make_exp_aluop_all_sizes(addm, (rd = ((ra + rb) & rc)), ZAPS);
make_exp_aluop_all_sizes(subm, (rd = ((ra - rb) & rc)), ZAPS);

make_exp_aluop_all_sizes(bt, (rb = lowbits(rb, log2(sizeof(T) * 8)), cf = bit(ra, rb), rd = (cf) ? -1 : +1), CF);
make_exp_aluop_all_sizes(bts, (rb = lowbits(rb, log2(sizeof(T) * 8)), cf = bit(ra, rb), rd = ra | (1LL << rb)), CF);
make_exp_aluop_all_sizes(btr, (rb = lowbits(rb, log2(sizeof(T) * 8)), cf = bit(ra, rb), rd = ra & ~(1LL << rb)), CF);
make_exp_aluop_all_sizes(btc, (rb = lowbits(rb, log2(sizeof(T) * 8)), cf = bit(ra, rb), rd = ra ^ (1LL << rb)), CF);

make_exp_aluop_all_sizes(bswap, (rd = ((sizeof(T) >= 4) ? std::byteswap(rb) : 0)), 0);

template<int ptlopcode, template<typename, int> class func, typename T, int genflags>
inline UopResult ctzclzop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  byte cf = 0, of = 0;
  func<T, genflags> f;
  T rt = f(ra, rb, rc, raflags, rbflags, rcflags, cf, of);
  out.rddata = x86_merge<T>(ra, rt);
  out.rdflags = (((T)rb) == 0) ? FLAG_ZF : 0;
  return out;
}

make_exp_aluop(exp_op_ctz, (rd = (rb) ? lsbindex64(rb) : 0));
make_anyop_all_sizes(OP_ctz, implmap_ctz, ctzclzop, exp_op_ctz, ZAPS);

make_exp_aluop(exp_op_clz, (rd = (rb) ? msbindex64(rb) : 0));
make_anyop_all_sizes(OP_clz, implmap_clz, ctzclzop, exp_op_clz, ZAPS);

UopResult uop_impl_collcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int flags = (raflags & FLAG_ZAPS) | (rbflags & FLAG_CF) | (rcflags & FLAG_OF);
  out.rddata = flags;
  out.rdflags = flags;
  return out;
}

UopResult uop_impl_movrcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int flags = rb & FLAG_NOT_WAIT_INV;
  out.rddata = flags;
  out.rdflags = flags;
  return out;
}

UopResult uop_impl_movccr(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int flags = rbflags;
  out.rddata = flags;
  out.rdflags = flags;
  return out;
}

UopResult uop_impl_andcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = (raflags & rbflags) & FLAG_NOT_WAIT_INV;
  return out;
}

UopResult uop_impl_orcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = (raflags | rbflags) & FLAG_NOT_WAIT_INV;
  return out;
}

UopResult uop_impl_ornotcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = (raflags | (~rbflags)) & FLAG_NOT_WAIT_INV;
  return out;
}

UopResult uop_impl_xorcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = (raflags ^ rbflags) & FLAG_NOT_WAIT_INV;
  return out;
}

make_x86_aluop_all_sizes(add, adc, ZAPS | CF | OF, PRETEXT_ALL_FLAGS_IN);
make_x86_aluop_all_sizes(sub, sbb, ZAPS | CF | OF, PRETEXT_ALL_FLAGS_IN);

//
// 3-operand ALU operation with shift of rc by 0/1/2/3
//

template<int ptlopcode, template<typename, int> class func, typename T, int genflags, int rcshift>
inline UopResult aluop3s(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  byte cf = 0, of = 0;
  func<T, genflags> f;
  T rt = f(ra, rb, rc << rcshift, raflags, rbflags, rcflags, cf, of);
  out.rddata = x86_merge<T>(ra, rt);
  // Do not generate of or cf for the 3-ops:
  out.rdflags = x86_genflags<T>(rt);
  return out;
}

// [size][extshift][setflags]
#define make_aluop3s_all_sizes_all_shifts(ptlopcode, mapname, nativeop, flagset)                                       \
  UopImpl mapname[4][4][2] = {                                                                                         \
      {                                                                                                                \
          {&aluop3s<ptlopcode, nativeop, W8, 0, 0>, &aluop3s<ptlopcode, nativeop, W8, (flagset), 0>},                  \
          {&aluop3s<ptlopcode, nativeop, W8, 0, 1>, &aluop3s<ptlopcode, nativeop, W8, (flagset), 1>},                  \
          {&aluop3s<ptlopcode, nativeop, W8, 0, 2>, &aluop3s<ptlopcode, nativeop, W8, (flagset), 2>},                  \
          {&aluop3s<ptlopcode, nativeop, W8, 0, 3>, &aluop3s<ptlopcode, nativeop, W8, (flagset), 3>},                  \
      },                                                                                                               \
      {                                                                                                                \
          {&aluop3s<ptlopcode, nativeop, W16, 0, 0>, &aluop3s<ptlopcode, nativeop, W16, (flagset), 0>},                \
          {&aluop3s<ptlopcode, nativeop, W16, 0, 1>, &aluop3s<ptlopcode, nativeop, W16, (flagset), 1>},                \
          {&aluop3s<ptlopcode, nativeop, W16, 0, 2>, &aluop3s<ptlopcode, nativeop, W16, (flagset), 2>},                \
          {&aluop3s<ptlopcode, nativeop, W16, 0, 3>, &aluop3s<ptlopcode, nativeop, W16, (flagset), 3>},                \
      },                                                                                                               \
      {                                                                                                                \
          {&aluop3s<ptlopcode, nativeop, W32, 0, 0>, &aluop3s<ptlopcode, nativeop, W32, (flagset), 0>},                \
          {&aluop3s<ptlopcode, nativeop, W32, 0, 1>, &aluop3s<ptlopcode, nativeop, W32, (flagset), 1>},                \
          {&aluop3s<ptlopcode, nativeop, W32, 0, 2>, &aluop3s<ptlopcode, nativeop, W32, (flagset), 2>},                \
          {&aluop3s<ptlopcode, nativeop, W32, 0, 3>, &aluop3s<ptlopcode, nativeop, W32, (flagset), 3>},                \
      },                                                                                                               \
      {                                                                                                                \
          {&aluop3s<ptlopcode, nativeop, W64, 0, 0>, &aluop3s<ptlopcode, nativeop, W64, (flagset), 0>},                \
          {&aluop3s<ptlopcode, nativeop, W64, 0, 1>, &aluop3s<ptlopcode, nativeop, W64, (flagset), 1>},                \
          {&aluop3s<ptlopcode, nativeop, W64, 0, 2>, &aluop3s<ptlopcode, nativeop, W64, (flagset), 2>},                \
          {&aluop3s<ptlopcode, nativeop, W64, 0, 3>, &aluop3s<ptlopcode, nativeop, W64, (flagset), 3>},                \
      },                                                                                                               \
  }

#define make_exp_aluop3_all_sizes_all_shifts(ptlopcode, name, expr, setflags)                                          \
  make_exp_aluop(exp_op_##name, (expr));                                                                               \
  make_aluop3s_all_sizes_all_shifts(ptlopcode, implmap_##name, exp_op_##name, (setflags));

make_exp_aluop3_all_sizes_all_shifts(OP_adda, adda, (rd = (ra + rb + rc)), 0);
make_exp_aluop3_all_sizes_all_shifts(OP_suba, suba, (rd = (ra - rb + rc)), 0);

//
// Shifts and rotates
//

template<typename T>
static inline unsigned int x86_shift_count(T rb) {
  const unsigned int width = sizeof(T) * 8;
  return byte(rb) & ((width == 64) ? 0x3f : 0x1f);
}

static inline byte flag_bit(W16 flags, W16 flag) {
  return (flags & flag) ? 1 : 0;
}

template<typename T>
static inline byte high_bit(T value) {
  using U = std::make_unsigned_t<T>;
  return byte((U(value) >> ((sizeof(T) * 8) - 1)) & 1);
}

template<typename U>
struct ShiftResult {
  U value;
  unsigned int count;
  byte cf;
  byte of;
};

struct ShiftLeftOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    const unsigned int count = x86_shift_count(rb);
    if (count == 0)
      return {value, count, 0, 0};

    const U result = (count < width) ? U(value << count) : U(0);
    const byte cf = (count <= width) ? byte((value >> (width - count)) & 1) : 0;
    const byte of = (count == 1) ? byte(high_bit(result) ^ cf) : 0;
    return {result, count, cf, of};
  }
};

struct ShiftRightOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    const unsigned int count = x86_shift_count(rb);
    if (count == 0)
      return {value, count, 0, 0};

    const U result = (count < width) ? U(value >> count) : U(0);
    const byte cf = (count <= width) ? byte((value >> (count - 1)) & 1) : 0;
    const byte of = (count == 1) ? high_bit(value) : 0;
    return {result, count, cf, of};
  }
};

struct ShiftArithmeticRightOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    const unsigned int count = x86_shift_count(rb);
    if (count == 0)
      return {value, count, 0, 0};

    const bool sign = high_bit(value) != 0;
    U result;
    if (count >= width)
      result = sign ? ~U(0) : U(0);
    else if (sign)
      result = U((value >> count) | (~U(0) << (width - count)));
    else
      result = U(value >> count);

    const byte cf = (count < width) ? byte((value >> (count - 1)) & 1) : byte(sign);
    return {result, count, cf, 0};
  }
};

struct RotateLeftOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    const unsigned int count = x86_shift_count(rb) % width;
    if (count == 0)
      return {value, count, 0, 0};

    const U result = std::rotl(value, int(count));
    const byte cf = byte(result & 1);
    const byte of = (count == 1) ? byte(high_bit(result) ^ cf) : 0;
    return {result, count, cf, of};
  }
};

struct RotateRightOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    const unsigned int count = x86_shift_count(rb) % width;
    if (count == 0)
      return {value, count, 0, 0};

    const U result = std::rotr(value, int(count));
    const byte cf = high_bit(result);
    const byte of = (count == 1) ? byte(((result >> (width - 1)) ^ (result >> (width - 2))) & 1) : 0;
    return {result, count, cf, of};
  }
};

struct RotateCarryLeftOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    unsigned int count = x86_shift_count(rb);
    if constexpr (sizeof(U) < sizeof(W32))
      count %= width + 1;
    if (count == 0)
      return {value, count, 0, 0};

    U result = value;
    byte cf = carry;
    for (unsigned int i = 0; i < count; i++) {
      byte nextcarry = high_bit(result);
      result = U((result << 1) | cf);
      cf = nextcarry;
    }

    const byte of = (count == 1) ? byte(high_bit(result) ^ cf) : 0;
    return {result, count, cf, of};
  }
};

struct RotateCarryRightOp {
  template<typename U, typename Count>
  ShiftResult<U> operator()(U value, Count rb, byte carry) const {
    const unsigned int width = sizeof(U) * 8;
    unsigned int count = x86_shift_count(rb);
    if constexpr (sizeof(U) < sizeof(W32))
      count %= width + 1;
    if (count == 0)
      return {value, count, 0, 0};

    U result = value;
    byte cf = carry;
    for (unsigned int i = 0; i < count; i++) {
      byte nextcarry = byte(result & 1);
      result = U((result >> 1) | (U(cf) << (width - 1)));
      cf = nextcarry;
    }

    const byte of = (count == 1) ? byte(((result >> (width - 1)) ^ (result >> (width - 2))) & 1) : 0;
    return {result, count, cf, of};
  }
};

template<int ptlopcode, typename Operation, typename T, int genflags>
inline UopResult shiftop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  using U = std::make_unsigned_t<T>;
  const auto result = Operation{}(U(ra), T(rb), flag_bit(rcflags, FLAG_CF));
  const T rt = T(result.value);
  out.rddata = x86_merge<T>(ra, rt);
  const byte cf = (genflags & SETFLAG_CF) ? result.cf : 0;
  const byte of = (genflags & SETFLAG_OF) ? result.of : 0;
  const int allflags = (of << 11) | cf | x86_genflags<T>(rt);
  out.rdflags = (result.count == 0) ? rcflags : allflags;
  return out;
}

template<int ptlopcode, typename Operation>
static constexpr std::array<std::array<UopImpl, 2>, 4> shiftop_impls = {
    std::array<UopImpl, 2>{&shiftop<ptlopcode, Operation, W8, 0>, &shiftop<ptlopcode, Operation, W8, ZAPS | CF | OF>},
    std::array<UopImpl, 2>{&shiftop<ptlopcode, Operation, W16, 0>, &shiftop<ptlopcode, Operation, W16, ZAPS | CF | OF>},
    std::array<UopImpl, 2>{&shiftop<ptlopcode, Operation, W32, 0>, &shiftop<ptlopcode, Operation, W32, ZAPS | CF | OF>},
    std::array<UopImpl, 2>{&shiftop<ptlopcode, Operation, W64, 0>,
                           &shiftop<ptlopcode, Operation, W64, ZAPS | CF | OF>}};

static constexpr auto implmap_shl = shiftop_impls<OP_shl, ShiftLeftOp>;
static constexpr auto implmap_shr = shiftop_impls<OP_shr, ShiftRightOp>;
static constexpr auto implmap_sar = shiftop_impls<OP_sar, ShiftArithmeticRightOp>;

static constexpr auto implmap_shls = shiftop_impls<OP_shls, ShiftLeftOp>;
static constexpr auto implmap_shrs = shiftop_impls<OP_shrs, ShiftRightOp>;
static constexpr auto implmap_sars = shiftop_impls<OP_sars, ShiftArithmeticRightOp>;

static constexpr auto implmap_rotl = shiftop_impls<OP_rotl, RotateLeftOp>;
static constexpr auto implmap_rotr = shiftop_impls<OP_rotr, RotateRightOp>;
static constexpr auto implmap_rotcl = shiftop_impls<OP_rotcl, RotateCarryLeftOp>;
static constexpr auto implmap_rotcr = shiftop_impls<OP_rotcr, RotateCarryRightOp>;

//
// Masks
//

template<int ptlopcode, typename T, int ZEROEXT, int SIGNEXT>
UopResult exp_op_mask(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  static const int sizeshift = log2(sizeof(T));
  W64 shmask = bitmask(6); //bitmask(3 + sizeshift);
  shmask = shmask | (shmask << 6) | (shmask << 12);
  rc &= shmask;

  int ms = bits(rc, 0, 6);
  int mc = bits(rc, 6, 6);
  int ds = bits(rc, 12, 6);

  int mcms = bits(rc, 0, 12);

  // M = 1'[(ms+mc-1):ms]
  W64 M = x86_ror<T>(static_cast<T>(bitmask(mc)), ms);
  W64 rd = (ra & ~M) | (x86_ror<T>(static_cast<T>(rb), ds) & M);

  logging::println(logging::TRACE,
                   "mask [{}, {}, {}, ss = {}, mcms {} [shmask {} (ms={} mc={} ds={} (mcms {}))]]:", sizeof(T), ZEROEXT,
                   SIGNEXT, sizeshift, mcms, bitstring(shmask, 18), ms, mc, ds, mcms);
  logging::println(logging::TRACE, "  M      = {} 0x{:016x}", bitstring(M, 64), M);
  logging::println(logging::TRACE, "  rot rb = {} 0x{:016x}", bitstring(x86_ror<T>(static_cast<T>(rb), ds), 64),
                   x86_ror<T>(static_cast<T>(rb), ds));
  logging::println(logging::TRACE, "  ra     = {:016x}", ra);
  logging::println(logging::TRACE, "  rb     = {:016x}", rb);
  logging::println(logging::TRACE, "  rc     = {:016x}", rc);
  logging::println(logging::TRACE, "  initrd = {:016x}", rd);

  if (ZEROEXT) {
    // rd = rd & 1'[(ms+mc-1):0]
    rd = rd & bitmask(ms + mc);
  } else if (SIGNEXT) {
    // rd = (rd[mc+ms-1]) ? (rd | 1'[63:(ms+mc)]) : (rd & 1'[(ms+mc-1):0]);
    rd = signext64(rd, ms + mc);
  } else {
    rd = rd;
  }

  out.rddata = x86_merge<T>(ra, rd);
  out.rdflags = x86_genflags<T>(rd);
  bool sf = bit(out.rdflags, log2(FLAG_SF));
  //
  // To simplify the microcode construction of the shrd instruction,
  // the following sequence may be used:
  //
  // shrd rd,rs:
  // shr  t = rd,c
  //      t.cf = rd[c-1] last bit shifted out
  //      t.of = rd[63]  or whatever rd's original sign bit position was
  // mask rd = t,rs,[ms=c, mc=c, ds=c]
  //      rd.cf = t.cf  inherited from t
  //      rd.of = (out.sf != t.of) i.e. did the sign bit change?
  //
  out.rdflags |= bit(raflags, log2(FLAG_CF)) << (log2(FLAG_CF));
  out.rdflags |= (sf != bit(raflags, log2(FLAG_OF))) << (log2(FLAG_OF));

  return out;
}

// [size][exttype]
UopImpl implmap_mask[4][3] = {
    {&exp_op_mask<OP_mask, W8, 0, 0>, &exp_op_mask<OP_mask, W8, 1, 0>, &exp_op_mask<OP_mask, W8, 0, 1>},
    {&exp_op_mask<OP_mask, W16, 0, 0>, &exp_op_mask<OP_mask, W16, 1, 0>, &exp_op_mask<OP_mask, W16, 0, 1>},
    {&exp_op_mask<OP_mask, W32, 0, 0>, &exp_op_mask<OP_mask, W32, 1, 0>, &exp_op_mask<OP_mask, W32, 0, 1>},
    {&exp_op_mask<OP_mask, W64, 0, 0>, &exp_op_mask<OP_mask, W64, 1, 0>, &exp_op_mask<OP_mask, W64, 0, 1>}};

// [size][exttype]
UopImpl implmap_maskb[4][3] = {
    {&exp_op_mask<OP_maskb, W8, 0, 0>, &exp_op_mask<OP_maskb, W8, 1, 0>, &exp_op_mask<OP_maskb, W8, 0, 1>},
    {&exp_op_mask<OP_maskb, W16, 0, 0>, &exp_op_mask<OP_maskb, W16, 1, 0>, &exp_op_mask<OP_maskb, W16, 0, 1>},
    {&exp_op_mask<OP_maskb, W32, 0, 0>, &exp_op_mask<OP_maskb, W32, 1, 0>, &exp_op_mask<OP_maskb, W32, 0, 1>},
    {&exp_op_mask<OP_maskb, W64, 0, 0>, &exp_op_mask<OP_maskb, W64, 1, 0>, &exp_op_mask<OP_maskb, W64, 0, 1>}};

//
// Permute bytes
//
// Technically this is a generalization of maskb, and maskb can be transformed
// into permb in the pipeline, at the cost of additional muxing logic.
//
UopResult uop_impl_permb(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  union vec128 {
    struct {
      W64 lo, hi;
    } w64;
    struct {
      byte b[16];
    } bytes;
  };

  union vec64 {
    struct {
      W64 data;
    } w64;
    struct {
      byte b[8];
    } bytes;
  };

  vec128 ab;
  vec64 d;

  ab.w64.lo = ra;
  ab.w64.hi = rb;

  logging::println(logging::TRACE, "Permute: by control 0x{:08x}:", rc);
  foreach (i, 8) {
    int which = bits(rc, i * 4, 4);

    logging::println(logging::TRACE, "  z[{}] = ab[{}] = 0x{:02x}", i, which, ab.bytes.b[which]);
    d.bytes.b[i] = ab.bytes.b[which];
  }

  out.rddata = d.w64.data;
  out.rdflags = x86_genflags<W64>(d.w64.data);

  return out;
}

//
// Multiplies
//

template<typename T>
struct double_width;
template<>
struct double_width<W8> {
  using u = W16;
  using s = W16s;
};
template<>
struct double_width<W16> {
  using u = W32;
  using s = W32s;
};
template<>
struct double_width<W32> {
  using u = W64;
  using s = W64s;
};
template<>
struct double_width<W64> {
  using u = unsigned __int128;
  using s = __int128;
};

enum class MulKind { Low, HighSigned, HighUnsigned };

template<MulKind kind, typename T, int genflags>
struct x86_mul {
  T operator()(T ra, T rb, T rc, W16 raflags, W16 rbflags, W16 rcflags, byte& cf, byte& of) {
    using S = std::make_signed_t<T>;
    using U2 = typename double_width<T>::u;
    using S2 = typename double_width<T>::s;
    constexpr int width = std::numeric_limits<T>::digits;

    const S2 sprod = S2(S(ra)) * S2(S(rb));
    const U2 uprod = U2(ra) * U2(rb);
    // One-operand MUL sets CF=OF if the high half is nonzero; one-operand IMUL
    // sets CF=OF if the full product is not the sign extension of the low half.
    cf = of = (kind == MulKind::HighUnsigned) ? byte((uprod >> width) != 0) : byte(sprod != S2(S(T(sprod))));
    if (kind == MulKind::Low)
      return T(sprod);
    return (kind == MulKind::HighSigned) ? T(U2(sprod) >> width) : T(uprod >> width);
  }
};

template<typename T, int genflags>
struct x86_op_mull : x86_mul<MulKind::Low, T, genflags> {};
template<typename T, int genflags>
struct x86_op_mulh : x86_mul<MulKind::HighSigned, T, genflags> {};
template<typename T, int genflags>
struct x86_op_mulhu : x86_mul<MulKind::HighUnsigned, T, genflags> {};

make_aluop_all_sizes(OP_mull, implmap_mull, x86_op_mull, ZAPS | CF | OF);
make_aluop_all_sizes(OP_mulh, implmap_mulh, x86_op_mulh, ZAPS | CF | OF);
make_aluop_all_sizes(OP_mulhu, implmap_mulhu, x86_op_mulhu, ZAPS | CF | OF);

// Can't fit 128-bit result in 64-bit register (otherwise it's the same as mull)
template<typename T, int genflags>
struct x86_op_mulhl {
  W64 operator()(W64 ra, W64 rb, W64 rc, W16 raflags, W16 rbflags, W16 rcflags, byte& cf, byte& of) {
    cf = 0;
    of = 0;
    T a = T(ra);
    T b = T(rb);
    return a * b;
  }
};

template<typename T>
inline UopResult uop_impl_mulhl(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T a = T(ra);
  T b = T(rb);
  W64 z = W64(a) * W64(b);

  /*
  // Do not merge: uop.size field then has ambiguous semantics
  W16 f = 0;
  switch (sizeof(T)) {
  case 1: z = x86_merge<W16>(ra, z); f = x86_genflags<W16>(z); break;
  case 2: z = x86_merge<W32>(ra, z); f = x86_genflags<W32>(z); break;
  case 4: z = x86_merge<W16>(ra, z); f = x86_genflags<W64>(z); break;
  case 8: z = z; f = x86_genflags<W64>(z); break;
  }
  */

  out.rddata = z;
  out.rdflags = x86_genflags<W64>(z);
  return out;
}

UopImpl implmap_mulhl[4] = {&uop_impl_mulhl<W8>, &uop_impl_mulhl<W16>, &uop_impl_mulhl<W32>, &uop_impl_mulhl<W64>};

template<int ptlopcode, typename T>
UopResult x86_div(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T quotient, remainder;

  if unlikely (!div_rem<T>(quotient, remainder, T(ra), T(rb), T(rc))) {
    out.rddata = EXCEPTION_DivideOverflow;
    out.rdflags = FLAG_INV;
    return out;
  }

  out.rddata = x86_merge<T>(rb, quotient);
  out.rdflags = x86_genflags<T>(quotient);
  return out;
}

template<int ptlopcode, typename T>
UopResult x86_rem(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T quotient, remainder;

  if unlikely (!div_rem<T>(quotient, remainder, T(ra), T(rb), T(rc))) {
    out.rddata = EXCEPTION_DivideOverflow;
    out.rdflags = FLAG_INV;
    return out;
  }

  out.rddata = x86_merge<T>(ra, remainder);
  out.rdflags = x86_genflags<T>(remainder);
  return out;
}

template<int ptlopcode, typename T>
UopResult x86_divs(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T quotient, remainder;

  if unlikely (!div_rem_s<T>(quotient, remainder, T(ra), T(rb), T(rc))) {
    out.rddata = EXCEPTION_DivideOverflow;
    out.rdflags = FLAG_INV;
    return out;
  }

  out.rddata = x86_merge<T>(rb, quotient);
  out.rdflags = x86_genflags<T>(quotient);
  return out;
}

template<int ptlopcode, typename T>
UopResult x86_rems(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T quotient, remainder;

  if unlikely (!div_rem_s<T>(quotient, remainder, T(ra), T(rb), T(rc))) {
    out.rddata = EXCEPTION_DivideOverflow;
    out.rdflags = FLAG_INV;
    return out;
  }

  out.rddata = x86_merge<T>(ra, remainder);
  out.rdflags = x86_genflags<T>(remainder);
  return out;
}

UopImpl implmap_div[4] = {&x86_div<OP_div, W8>, &x86_div<OP_div, W16>, &x86_div<OP_div, W32>, &x86_div<OP_div, W64>};
UopImpl implmap_rem[4] = {&x86_rem<OP_rem, W8>, &x86_rem<OP_rem, W16>, &x86_rem<OP_rem, W32>, &x86_rem<OP_rem, W64>};
UopImpl implmap_divs[4] = {&x86_divs<OP_divs, W8>, &x86_divs<OP_divs, W16>, &x86_divs<OP_divs, W32>,
                           &x86_divs<OP_divs, W64>};
UopImpl implmap_rems[4] = {&x86_rems<OP_rems, W8>, &x86_rems<OP_rems, W16>, &x86_rems<OP_rems, W32>,
                           &x86_rems<OP_rems, W64>};

template<int ptlopcode, typename T, bool compare_for_max>
UopResult uop_impl_min_max(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  T a = ra;
  T b = rb;
  T z = (compare_for_max) ? std::max(a, b) : std::min(a, b);
  out.rddata = x86_merge<T>(ra, z);
  out.rdflags = x86_genflags<T>(z);
  return out;
}

UopImpl implmap_min[4] = {&uop_impl_min_max<OP_min, W8, 0>, &uop_impl_min_max<OP_min, W16, 0>,
                          &uop_impl_min_max<OP_min, W32, 0>, &uop_impl_min_max<OP_min, W64, 0>};
UopImpl implmap_max[4] = {&uop_impl_min_max<OP_max, W8, 1>, &uop_impl_min_max<OP_max, W16, 1>,
                          &uop_impl_min_max<OP_max, W32, 1>, &uop_impl_min_max<OP_max, W64, 1>};
UopImpl implmap_min_s[4] = {&uop_impl_min_max<OP_min, W8s, 0>, &uop_impl_min_max<OP_min, W16s, 0>,
                            &uop_impl_min_max<OP_min, W32s, 0>, &uop_impl_min_max<OP_min, W64s, 0>};
UopImpl implmap_max_s[4] = {&uop_impl_min_max<OP_max, W8s, 1>, &uop_impl_min_max<OP_max, W16s, 1>,
                            &uop_impl_min_max<OP_max, W32s, 1>, &uop_impl_min_max<OP_max, W64s, 1>};

//
// Condition code evaluation
//

template<int evaltype>
inline bool evaluate_cond(int ra, int rb) {
  switch (evaltype) {
  case 0: // {0, REG_zero, REG_of},   // of:               jo
    return !!(rb & FLAG_OF);
  case 1: // {0, REG_zero, REG_of},   // !of:              jno
    return !(rb & FLAG_OF);
  case 2: // {0, REG_zero, REG_cf},   // cf:               jb jc jnae
    return !!(rb & FLAG_CF);
  case 3: // {0, REG_zero, REG_cf},   // !cf:              jnb jnc jae
    return !(rb & FLAG_CF);
  case 4: // {0, REG_zf,   REG_zero}, // zf:               jz je
    return !!(ra & FLAG_ZF);
  case 5: // {0, REG_zf,   REG_zero}, // !zf:              jnz jne
    return !(ra & FLAG_ZF);
  case 6: // {1, REG_zf,   REG_cf},   // cf|zf:            jbe jna
    return ((ra & FLAG_ZF) || (rb & FLAG_CF));
  case 7: // {1, REG_zf,   REG_cf},   // !cf & !zf:        jnbe ja
    return !((ra & FLAG_ZF) || (rb & FLAG_CF));
  case 8: // {0, REG_zf,   REG_zero}, // sf:               js
    return !!(ra & FLAG_SF);
  case 9: // {0, REG_zf,   REG_zero}, // !sf:              jns
    return !(ra & FLAG_SF);
  case 10: // {0, REG_zf,   REG_zero}, // pf:               jp jpe
    return !!(ra & FLAG_PF);
  case 11: // {0, REG_zf,   REG_zero}, // !pf:              jnp jpo
    return !(ra & FLAG_PF);
  case 12: // {1, REG_zf,   REG_of},   // sf != of:         jl jnge (*)
    return (!!(ra & FLAG_SF)) != (!!(rb & FLAG_OF));
  case 13: // {1, REG_zf,   REG_of},   // sf == of:         jnl jge (*)
    return !(!!(ra & FLAG_SF)) != (!!(rb & FLAG_OF));
  case 14: // {1, REG_zf,   REG_of},   // zf | (sf != of):  jle jng (*)
    return ((!!(ra & FLAG_ZF)) | ((!!(ra & FLAG_SF)) != (!!(rb & FLAG_OF))));
  case 15: // {1, REG_zf,   REG_of},   // !zf & (sf == of): jnle jg (*)
    return !((!!(ra & FLAG_ZF)) | ((!!(ra & FLAG_SF)) != (!!(rb & FLAG_OF))));
  }
}

#define make_condop_all_conds_any(ptlopcode, subtype, subarrays, mapname, operation)                                   \
  UopImpl implmap_##mapname[16] subarrays = {                                                                          \
      subtype(ptlopcode, operation, 0),  subtype(ptlopcode, operation, 1),  subtype(ptlopcode, operation, 2),          \
      subtype(ptlopcode, operation, 3),  subtype(ptlopcode, operation, 4),  subtype(ptlopcode, operation, 5),          \
      subtype(ptlopcode, operation, 6),  subtype(ptlopcode, operation, 7),  subtype(ptlopcode, operation, 8),          \
      subtype(ptlopcode, operation, 9),  subtype(ptlopcode, operation, 10), subtype(ptlopcode, operation, 11),         \
      subtype(ptlopcode, operation, 12), subtype(ptlopcode, operation, 13), subtype(ptlopcode, operation, 14),         \
      subtype(ptlopcode, operation, 15)}

#define make_condop(ptlopcode, operation, cond) &operation<ptlopcode, cond>
#define make_condop_all_sizes(ptlopcode, operation, cond)                                                              \
  {&operation<ptlopcode, W8, cond>, &operation<ptlopcode, W16, cond>, &operation<ptlopcode, W32, cond>,                \
   &operation<ptlopcode, W64, cond>}

#define make_condop_all_conds_all_sizes(mapname, operation)                                                            \
  make_condop_all_conds_any(OP_##mapname, make_condop_all_sizes, [4], mapname, operation)

#define function(expr, rettype, ...)                                                                                   \
  class {                                                                                                              \
  public:                                                                                                              \
    rettype operator()(__VA_ARGS__) { return (expr); }                                                                 \
  }

template<typename T>
struct sub_flag_gen_op {
  W16 operator()(T ra, T rb) {
    x86_op_sub<T, ZAPS | CF | OF> op;
    byte cf, of;
    T rd = op(ra, rb, 0, 0, 0, 0, cf, of);
    return (of << 11) | cf | x86_genflags<T>(rd);
  }
};

template<typename T>
struct and_flag_gen_op {
  W16 operator()(T ra, T rb) { return x86_genflags<T>(ra & rb); }
};

//
// sel.cc, sel.cmp.cc
//
template<int ptlopcode, typename T, int evaltype>
inline UopResult uop_impl_sel(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool istrue = evaluate_cond<evaltype>(rcflags, rcflags);
  out.rddata = x86_merge<T>(ra, (istrue) ? rb : ra);
  out.rdflags = (istrue) ? rbflags : raflags;
  return out;
}

make_condop_all_conds_all_sizes(sel, uop_impl_sel);

#define make_condop_all_sizes_all_compare_sizes(ptlopcode, operation, cond)                                            \
  {                                                                                                                    \
      {&operation<ptlopcode, W8, W8, cond>, &operation<ptlopcode, W8, W16, cond>,                                      \
       &operation<ptlopcode, W8, W32, cond>, &operation<ptlopcode, W8, W64, cond>},                                    \
      {&operation<ptlopcode, W16, W8, cond>, &operation<ptlopcode, W16, W16, cond>,                                    \
       &operation<ptlopcode, W16, W32, cond>, &operation<ptlopcode, W16, W64, cond>},                                  \
      {&operation<ptlopcode, W32, W8, cond>, &operation<ptlopcode, W32, W16, cond>,                                    \
       &operation<ptlopcode, W32, W32, cond>, &operation<ptlopcode, W32, W64, cond>},                                  \
      {&operation<ptlopcode, W64, W8, cond>, &operation<ptlopcode, W64, W16, cond>,                                    \
       &operation<ptlopcode, W64, W32, cond>, &operation<ptlopcode, W64, W64, cond>},                                  \
  }

#define make_condop_all_conds_all_sizes_all_compare_sizes(ptlopcode, operation)                                        \
  make_condop_all_conds_any(OP_##ptlopcode, make_condop_all_sizes_all_compare_sizes, [4][4], ptlopcode, operation)

template<int ptlopcode, typename Tmerge, typename Tcompare, int evaltype>
inline UopResult uop_impl_sel_cmp(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int flags = x86_genflags<Tcompare>(rc);
  bool istrue = evaluate_cond<evaltype>(flags, flags);
  out.rddata = x86_merge<Tmerge>(ra, (istrue) ? rb : ra);
  out.rdflags = (istrue) ? rbflags : raflags;
  return out;
}

make_condop_all_conds_all_sizes_all_compare_sizes(sel_cmp, uop_impl_sel_cmp);

//
// set.cc, set.sub.cc, set.and.cc
//
template<int ptlopcode, typename T, int evaltype>
inline UopResult uop_impl_set(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool istrue = evaluate_cond<evaltype>(raflags, rbflags);
  out.rddata = x86_merge<T>(rc, (istrue) ? 1 : 0);
  out.rdflags = (istrue) ? FLAG_CF : 0;
  return out;
}

make_condop_all_conds_all_sizes(set, uop_impl_set);

template<int ptlopcode, typename Tmerge, typename Tcompare, int evaltype>
inline UopResult uop_impl_set_and(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  and_flag_gen_op<Tcompare> func;
  int flags = func(ra, rb);
  bool istrue = evaluate_cond<evaltype>(flags, flags);
  out.rddata = x86_merge<Tmerge>(rc, (istrue) ? 1 : 0);
  out.rdflags = (istrue) ? FLAG_CF : 0;
  return out;
}

make_condop_all_conds_all_sizes_all_compare_sizes(set_and, uop_impl_set_and);

template<int ptlopcode, typename Tmerge, typename Tcompare, int evaltype>
inline UopResult uop_impl_set_sub(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  sub_flag_gen_op<Tcompare> func;
  int flags = func(ra, rb);
  bool istrue = evaluate_cond<evaltype>(flags, flags);
  out.rddata = x86_merge<Tmerge>(rc, (istrue) ? 1 : 0);
  out.rdflags = (istrue) ? FLAG_CF : 0;
  return out;
}

make_condop_all_conds_all_sizes_all_compare_sizes(set_sub, uop_impl_set_sub);

//
// Branches
//

template<int ptlopcode, typename T, int evaltype, bool excepting>
inline UopResult uop_impl_condbranch(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool taken = evaluate_cond<evaltype>(raflags, rbflags);
  out.rddata = (taken) ? riptaken : ripseq;
  out.rdflags = (taken) ? FLAG_BR_TK : 0;

  if (excepting & (!taken)) {
    out.rddata = EXCEPTION_BranchMispredict;
    out.rdflags |= FLAG_INV;
  }
  return out;
}

#define make_branchop_all_excepts(ptlopcode, operation, cond)                                                          \
  {&uop_impl_condbranch<ptlopcode, W64, cond, false>, &uop_impl_condbranch<ptlopcode, W64, cond, true>}

make_condop_all_conds_any(OP_br, make_branchop_all_excepts, [2], br, anything);

template<int ptlopcode, typename T, int evaltype, bool excepting, template<typename> class func_t>
inline UopResult uop_impl_alu_and_condbranch(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  func_t<T> func;
  int flags = func(ra, rb);
  bool taken = evaluate_cond<evaltype>(flags, flags);
  out.rddata = (taken) ? riptaken : ripseq;
  out.rdflags = flags | (taken ? FLAG_BR_TK : 0);

  if (excepting & (!taken)) {
    out.rddata = EXCEPTION_BranchMispredict;
    out.rdflags |= FLAG_INV;
  }
  return out;
}

#define make_alu_and_branchop_all_sizes_all_excepts(ptlopcode, operation, cond)                                        \
  {                                                                                                                    \
      {&uop_impl_alu_and_condbranch<ptlopcode, W8, cond, false, operation>,                                            \
       &uop_impl_alu_and_condbranch<ptlopcode, W8, cond, true, operation>},                                            \
      {&uop_impl_alu_and_condbranch<ptlopcode, W16, cond, false, operation>,                                           \
       &uop_impl_alu_and_condbranch<ptlopcode, W16, cond, true, operation>},                                           \
      {&uop_impl_alu_and_condbranch<ptlopcode, W32, cond, false, operation>,                                           \
       &uop_impl_alu_and_condbranch<ptlopcode, W32, cond, true, operation>},                                           \
      {&uop_impl_alu_and_condbranch<ptlopcode, W64, cond, false, operation>,                                           \
       &uop_impl_alu_and_condbranch<ptlopcode, W64, cond, true, operation>},                                           \
  }

make_condop_all_conds_any(OP_br_and, make_alu_and_branchop_all_sizes_all_excepts, [4][2], br_and, and_flag_gen_op);
make_condop_all_conds_any(OP_br_sub, make_alu_and_branchop_all_sizes_all_excepts, [4][2], br_sub, sub_flag_gen_op);

UopResult uop_impl_jmp(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool taken = (riptaken == ra);
  out.rddata = ra;
  out.rdflags = (taken) ? FLAG_BR_TK : 0;
  return out;
}

UopResult uop_impl_jmp_ex(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool taken = (riptaken == ra);
  out.rddata = ra;
  out.rdflags = (taken) ? FLAG_BR_TK : 0;

  if (!taken) {
    out.rddata = EXCEPTION_BranchMispredict;
    out.rdflags |= FLAG_INV;
  }
  return out;
}

UopResult uop_impl_bru(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = riptaken;
  out.rdflags = FLAG_BR_TK;
  return out;
}

UopResult uop_impl_brp(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = riptaken;
  out.rdflags = FLAG_BR_TK;
  return out;
}

//
// Checks
//
template<int ptlopcode, int evaltype>
inline UopResult uop_impl_chk(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  bool passed = evaluate_cond<evaltype>(raflags, rbflags);
  out.rddata = (passed) ? 0 : rc;
  out.addr = 0;
  out.rdflags = (passed) ? 0 : FLAG_INV;
  return out;
}

template<int ptlopcode, typename T, int evaltype>
inline UopResult uop_impl_chk_sub(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  sub_flag_gen_op<T> func;
  int flags = func(ra, rb);
  bool passed = evaluate_cond<evaltype>(flags, flags);
  out.rddata = (passed) ? 0 : rc;
  out.addr = 0;
  out.rdflags = (passed) ? 0 : FLAG_INV;
  return out;
}

template<int ptlopcode, typename T, int evaltype>
inline UopResult uop_impl_chk_and(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  and_flag_gen_op<T> func;
  int flags = func(ra, rb);
  bool passed = evaluate_cond<evaltype>(flags, flags);
  out.rddata = (passed) ? 0 : rc;
  out.addr = 0;
  out.rdflags = (passed) ? 0 : FLAG_INV;
  return out;
}

make_condop_all_conds_any(OP_chk, make_condop, [1], chk, uop_impl_chk);
make_condop_all_conds_all_sizes(chk_sub, uop_impl_chk_sub);
make_condop_all_conds_all_sizes(chk_and, uop_impl_chk_and);


//
// Floating Point
//
#define make_exp_floatop(name, expr)                                                                                   \
  template<typename T>                                                                                                 \
  struct name {                                                                                                        \
    T operator()(T ra, T rb, T rc) {                                                                                   \
      T rd;                                                                                                            \
      expr;                                                                                                            \
      return rd;                                                                                                       \
    }                                                                                                                  \
  }

template<int ptlopcode, template<typename> class F, int datatype>
inline UopResult floatop(const UopInputs& inputs) {
  [[maybe_unused]] auto [raraw, rbraw, rcraw, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  SSEType ra, rb, rc, rd;
  ra.w64 = raraw;
  rb.w64 = rbraw;
  rc.w64 = rcraw;

  switch (datatype) {
  case 0: { // scalar single
    F<float> func;
    rd.f.lo = func(ra.f.lo, rb.f.lo, rc.f.lo);
    rd.w32.hi = ra.w32.hi;
    break;
  }
  case 1: { // packed single
    F<float> func;
    rd.f.lo = func(ra.f.lo, rb.f.lo, rc.f.lo);
    rd.f.hi = func(ra.f.hi, rb.f.hi, rc.f.hi);
    break;
  }
  case 2:
  case 3: { // scalar double
    F<float> func;
    rd.d = func(ra.d, rb.d, rc.d);
    break;
  }
  }
  out.rddata = rd.w64;
  out.rdflags = 0;
  return out;
}

#define make_exp_floatop_alltypes(name, expr)                                                                          \
  make_exp_floatop(exp_op_##name, expr);                                                                               \
  UopImpl implmap_##name[4] = {&floatop<OP_##name, exp_op_##name, 0>, &floatop<OP_##name, exp_op_##name, 1>,           \
                               &floatop<OP_##name, exp_op_##name, 2>, &floatop<OP_##name, exp_op_##name, 3>}

enum class SSEFloatType : int {
  ScalarSingle = 0,
  PackedSingle = 1,
  ScalarDouble = 2,
  PackedDouble = 3,
};

template<int predicate, typename T>
static inline T fcmp_result(T ra, T rb) {
  using Bits = std::conditional_t<sizeof(T) == sizeof(W32), W32, W64>;

  const bool unordered = std::isunordered(ra, rb);
  bool result;

  if constexpr (predicate == 0)
    result = (ra == rb);
  else if constexpr (predicate == 1)
    result = (ra < rb);
  else if constexpr (predicate == 2)
    result = (ra <= rb);
  else if constexpr (predicate == 3)
    result = unordered;
  else if constexpr (predicate == 4)
    result = (ra != rb);
  else if constexpr (predicate == 5)
    result = !(ra < rb);
  else if constexpr (predicate == 6)
    result = !(ra <= rb);
  else
    result = !unordered;

  T all_ones = std::bit_cast<T>(~Bits(0));
  return result ? all_ones : T(0);
}

template<int ptlopcode, SSEFloatType type, typename Operation>
UopResult floatop2(const UopInputs& inputs) {
  [[maybe_unused]] auto [raraw, rbraw, rcraw, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  SSEType ra, rb, rd;
  ra.w64 = raraw;
  rb.w64 = rbraw;

  Operation operation;

  if constexpr (type == SSEFloatType::ScalarSingle) {
    rd.f.lo = operation(ra.f.lo, rb.f.lo);
    rd.w32.hi = ra.w32.hi;
  } else if constexpr (type == SSEFloatType::PackedSingle) {
    rd.f.lo = operation(ra.f.lo, rb.f.lo);
    rd.f.hi = operation(ra.f.hi, rb.f.hi);
  } else {
    rd.d = operation(ra.d, rb.d);
  }

  out.rddata = rd.w64;
  out.rdflags = 0;
  return out;
}

struct FloatSqrt {
  template<typename T>
  T operator()(T, T rb) const {
    return std::sqrt(rb);
  }
};

struct FloatReciprocal {
  template<typename T>
  T operator()(T, T rb) const {
    return std::divides<>{}(T(1), rb);
  }
};

struct FloatReciprocalSqrt {
  template<typename T>
  T operator()(T, T rb) const {
    return std::divides<>{}(T(1), std::sqrt(rb));
  }
};

struct FloatMin {
  template<typename T>
  T operator()(T ra, T rb) const {
    if (std::isnan(ra) || std::isnan(rb) || (ra == rb))
      return rb;

    return std::min(ra, rb);
  }
};

struct FloatMax {
  template<typename T>
  T operator()(T ra, T rb) const {
    if (std::isnan(ra) || std::isnan(rb) || (ra == rb))
      return rb;

    return std::max(ra, rb);
  }
};

template<int predicate>
struct FloatCompare {
  template<typename T>
  T operator()(T ra, T rb) const {
    return fcmp_result<predicate>(ra, rb);
  }
};

template<int ptlopcode, typename Operation>
static constexpr std::array<UopImpl, 4> floatop2_impls = {&floatop2<ptlopcode, SSEFloatType::ScalarSingle, Operation>,
                                                          &floatop2<ptlopcode, SSEFloatType::PackedSingle, Operation>,
                                                          &floatop2<ptlopcode, SSEFloatType::ScalarDouble, Operation>,
                                                          &floatop2<ptlopcode, SSEFloatType::PackedDouble, Operation>};

static constexpr auto implmap_fadd = floatop2_impls<OP_fadd, std::plus<>>;
static constexpr auto implmap_fsub = floatop2_impls<OP_fsub, std::minus<>>;
static constexpr auto implmap_fmul = floatop2_impls<OP_fmul, std::multiplies<>>;
static constexpr auto implmap_fdiv = floatop2_impls<OP_fdiv, std::divides<>>;
static constexpr auto implmap_fsqrt = floatop2_impls<OP_fsqrt, FloatSqrt>;
static constexpr auto implmap_frcp = floatop2_impls<OP_frcp, FloatReciprocal>;
static constexpr auto implmap_frsqrt = floatop2_impls<OP_frsqrt, FloatReciprocalSqrt>;
static constexpr auto implmap_fmin = floatop2_impls<OP_fmin, FloatMin>;
static constexpr auto implmap_fmax = floatop2_impls<OP_fmax, FloatMax>;

//
// Derived operations:
//
// fadd   +(ra * 1.0) + rc     =>  fmadd ra, 1.0, rc
// fsub   +(ra * 1.0) - rc     =>  fmsub ra, 1.0, rc
// fmul   +(ra * rb)  - 0.0    =>  fmadd ra, rb,  -0.0
//

template<SSEFloatType type>
UopResult x86_op_fmadd(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  //
  // fmadd  rd = (ra * rb) + rc       =>  fmul t0 = ra,rb  |  fadd t1 = t0,rc
  //
  const UopResult product = floatop2<OP_fmul, type, std::multiplies<>>({.ra = ra, .rb = rb});
  return floatop2<OP_fsub, type, std::plus<>>({.ra = product.rddata, .rb = rc});
}

UopImpl implmap_fmadd[4] = {&x86_op_fmadd<SSEFloatType::ScalarSingle>, &x86_op_fmadd<SSEFloatType::PackedSingle>,
                            &x86_op_fmadd<SSEFloatType::ScalarDouble>, &x86_op_fmadd<SSEFloatType::PackedDouble>};

template<SSEFloatType type>
UopResult x86_op_fmsub(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  //
  // fmsub  rd = (ra * rb) - rc       =>  fmul t0 = ra,rb  |  fsub t1 = t0,rc
  //
  const UopResult product = floatop2<OP_fmul, type, std::multiplies<>>({.ra = ra, .rb = rb});
  return floatop2<OP_fsub, type, std::minus<>>({.ra = product.rddata, .rb = rc});
}

UopImpl implmap_fmsub[4] = {&x86_op_fmsub<SSEFloatType::ScalarSingle>, &x86_op_fmsub<SSEFloatType::PackedSingle>,
                            &x86_op_fmsub<SSEFloatType::ScalarDouble>, &x86_op_fmsub<SSEFloatType::PackedDouble>};

template<SSEFloatType type>
UopResult x86_op_fmsubr(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  //
  // fmsubr rd = rc - (ra * rb)       =>  fmul t0 = ra,rb  |  fsub t1 = rc,t0
  //
  const UopResult product = floatop2<OP_fmul, type, std::multiplies<>>({.ra = ra, .rb = rb});
  return floatop2<OP_fsub, type, std::minus<>>({.ra = rc, .rb = product.rddata});
}

UopImpl implmap_fmsubr[4] = {&x86_op_fmsubr<SSEFloatType::ScalarSingle>, &x86_op_fmsubr<SSEFloatType::PackedSingle>,
                             &x86_op_fmsubr<SSEFloatType::ScalarDouble>, &x86_op_fmsubr<SSEFloatType::PackedDouble>};

static constexpr std::array<std::array<UopImpl, 4>, 8> implmap_fcmp = {
    floatop2_impls<OP_fcmp, FloatCompare<0>>, floatop2_impls<OP_fcmp, FloatCompare<1>>,
    floatop2_impls<OP_fcmp, FloatCompare<2>>, floatop2_impls<OP_fcmp, FloatCompare<3>>,
    floatop2_impls<OP_fcmp, FloatCompare<4>>, floatop2_impls<OP_fcmp, FloatCompare<5>>,
    floatop2_impls<OP_fcmp, FloatCompare<6>>, floatop2_impls<OP_fcmp, FloatCompare<7>>};

// comis/ucomis (comptype 0/1: float, 2/3: double): ra is the x86 destination
// operand, rb the source. SIMD FP exceptions are not modeled, so the ordered
// and unordered variants behave identically.
template<int comptype>
UopResult uop_impl_fcmpcc(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  constexpr bool dbl = (comptype >= 2);
  const SSEType a(ra), b(rb);
  const double x = dbl ? a.d : a.f.lo;
  const double y = dbl ? b.d : b.f.lo;

  byte zf, pf, cf;
  if (std::isunordered(x, y)) {
    zf = pf = cf = 1;
  } else {
    zf = (x == y);
    pf = 0;
    cf = (x < y);
  }
  out.rdflags = (zf << 6) + (pf << 2) + (cf << 0);
  out.rddata = out.rdflags;
  return out;
}

UopImpl implmap_fcmpcc[8][4] = {&uop_impl_fcmpcc<0>, &uop_impl_fcmpcc<1>, &uop_impl_fcmpcc<2>, &uop_impl_fcmpcc<3>};

#define make_intsrc_fp_convop(name, op)                                                                                \
  UopResult uop_impl_##name(const UopInputs& inputs) {                                                                 \
    [[maybe_unused]] auto [raraw, rbraw, rcraw, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;                 \
    UopResult out;                                                                                                     \
    SSEType ra, rb, rc, rd;                                                                                            \
    ra.w64 = raraw;                                                                                                    \
    rb.w64 = rbraw;                                                                                                    \
    rc.w64 = rcraw;                                                                                                    \
    op;                                                                                                                \
    out.rddata = rd.w64;                                                                                               \
    out.rdflags = 0;                                                                                                   \
    return out;                                                                                                        \
  }

make_intsrc_fp_convop(fcvt_i2s_ins, (rd.f.lo = (float)(W32s)rb.w32.lo, rd.w32.hi = ra.w32.hi));
make_intsrc_fp_convop(fcvt_q2s_ins, (rd.f.lo = (float)(W64s)rb.w64, rd.w32.hi = ra.w32.hi));
make_intsrc_fp_convop(fcvt_q2d, (rd.d = (double)(W64s)rb.w64));

make_intsrc_fp_convop(fcvt_i2s_p, (rd.f.lo = (float)(W32s)rb.w32.lo, rd.f.hi = (float)(W32s)rb.w32.hi));
make_intsrc_fp_convop(fcvt_i2d_lo, (rd.d = (double)(W32s)rb.w32.lo));
make_intsrc_fp_convop(fcvt_i2d_hi, (rd.d = (double)(W32s)rb.w32.hi));
make_intsrc_fp_convop(fcvt_s2d_lo, (rd.d = (double)rb.f.lo));
make_intsrc_fp_convop(fcvt_s2d_hi, (rd.d = (double)rb.f.hi));
make_intsrc_fp_convop(fcvt_d2s_ins, (rd.f.lo = (float)rb.d, rd.w32.hi = ra.w32.hi));

//
// Convert a scalar FP value to a signed integer with x86 SSE cvt/cvtt
// semantics: the trunc variant rounds toward zero, the rounding variant uses
// round to nearest even (the default of the unmodeled MXCSR); NaN and
// out-of-range inputs yield the "integer indefinite" value (1 << (bits-1)).
// float -> double is exact, so one double-based helper covers both widths.
//
template<typename Int>
static inline Int x86_fp_to_int(double v, bool trunc) {
  const double r = trunc ? std::trunc(v) : std::nearbyint(v);
  constexpr double limit = (sizeof(Int) == 4) ? 0x1p31 : 0x1p63;
  if (!(r >= -limit && r < limit))
    return Int(1) << (std::numeric_limits<Int>::digits - 1);
  return Int(std::make_signed_t<Int>(r));
}

#define make_intdest_fp_convop_allrounds(name, desttype, srcexpr)                                                      \
  template<int ptlopcode, int trunc>                                                                                   \
  UopResult uop_impl_##name(const UopInputs& inputs) {                                                                 \
    [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;                          \
    UopResult out;                                                                                                     \
    const SSEType src(rb);                                                                                             \
    out.rddata = x86_fp_to_int<desttype>(srcexpr, trunc);                                                              \
    out.rdflags = 0;                                                                                                   \
    return out;                                                                                                        \
  }                                                                                                                    \
  UopImpl implmap_##name[2] = {&uop_impl_##name<OP_##name, 0>, &uop_impl_##name<OP_##name, 1>}

make_intdest_fp_convop_allrounds(fcvt_s2i, W32, src.f.lo);
make_intdest_fp_convop_allrounds(fcvt_d2i, W32, src.d);
make_intdest_fp_convop_allrounds(fcvt_s2q, W64, src.f.lo);
make_intdest_fp_convop_allrounds(fcvt_d2q, W64, src.d);

#define make_fp_convop_allrounds(name, expr)                                                                           \
  template<int trunc>                                                                                                  \
  UopResult uop_impl_##name(const UopInputs& inputs) {                                                                 \
    [[maybe_unused]] auto [raraw, rbraw, rcraw, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;                 \
    UopResult out;                                                                                                     \
    const SSEType ra(raraw), rb(rbraw);                                                                                \
    SSEType rd;                                                                                                        \
    expr;                                                                                                              \
    out.rddata = rd.w64;                                                                                               \
    out.rdflags = 0;                                                                                                   \
    return out;                                                                                                        \
  }                                                                                                                    \
  UopImpl implmap_##name[2] = {&uop_impl_##name<0>, &uop_impl_##name<1>}

// Note the lane routing (from the original cvt[t]pd2dq/cvtpd2ps operand setup):
// the packed-double ops take the low lane from rb and the high lane from ra.
make_fp_convop_allrounds(fcvt_s2i_p, (rd.w32.lo = x86_fp_to_int<W32>(rb.f.lo, trunc),
                                      rd.w32.hi = x86_fp_to_int<W32>(rb.f.hi, trunc)));
make_fp_convop_allrounds(fcvt_d2i_p,
                         (rd.w32.lo = x86_fp_to_int<W32>(rb.d, trunc), rd.w32.hi = x86_fp_to_int<W32>(ra.d, trunc)));
make_fp_convop_allrounds(fcvt_d2s_p, (rd.f.lo = (float)rb.d, rd.f.hi = (float)ra.d));

//
// Vector uops (MMX width: one 64-bit register holding packed lanes)
//

// Dummy (to fill in unsupported places only)
template<int ptlopcode, int sizeshift>
UopResult x86_op_nop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  out.rddata = 0;
  out.rdflags = 0;
  return out;
}

template<typename Lane>
using LaneArray = std::array<Lane, 8 / sizeof(Lane)>;

// Saturate a wide intermediate to the full range of Lane
template<typename Lane>
static constexpr Lane saturate(W64s v) {
  return Lane(std::clamp<W64s>(v, std::numeric_limits<Lane>::min(), std::numeric_limits<Lane>::max()));
}

template<int ptlopcode, typename Lane, auto op>
UopResult vecop(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  auto a = std::bit_cast<LaneArray<Lane>>(ra);
  const auto b = std::bit_cast<LaneArray<Lane>>(rb);
  foreach (i, a.size())
    a[i] = Lane(op(a[i], b[i]));
  out.rddata = std::bit_cast<W64>(a);
  out.rdflags = 0;
  return out;
}

// MMX shifts take the entire 64-bit rb as one shift count, deliberately
// unmasked: counts >= the lane width zero the lanes (or sign-fill for psra).
enum class VecShift { Left, Right, RightArith };

template<int ptlopcode, typename Lane, VecShift kind>
UopResult vecshift(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  constexpr W64 width = sizeof(Lane) * 8;
  auto a = std::bit_cast<LaneArray<Lane>>(ra);
  foreach (i, a.size()) {
    if constexpr (kind == VecShift::Left)
      a[i] = (rb < width) ? Lane(a[i] << rb) : Lane(0);
    else if constexpr (kind == VecShift::Right)
      a[i] = (rb < width) ? Lane(a[i] >> rb) : Lane(0);
    else
      a[i] = Lane(std::make_signed_t<Lane>(a[i]) >> std::min(rb, width - 1));
  }
  out.rddata = std::bit_cast<W64>(a);
  out.rdflags = 0;
  return out;
}

// Narrow 2N source lanes to N-wide destination lanes with saturation; the low
// half of the result comes from ra, the high half from rb.
template<int ptlopcode, typename Lane, bool unsigned_sat>
UopResult vecpack(const UopInputs& inputs) {
  [[maybe_unused]] auto [raraw, rbraw, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  // The guard keeps the wide (unused, nop-dispatched) lane widths from being
  // instantiated by make_vec_implmap: their source lanes would exceed 64 bits.
  if constexpr (sizeof(Lane) <= 2) {
    using Src = typename double_width<Lane>::s;
    using SLane = std::make_signed_t<Lane>;
    const auto a = std::bit_cast<LaneArray<Src>>(raraw);
    const auto b = std::bit_cast<LaneArray<Src>>(rbraw);
    LaneArray<Lane> d;
    foreach (i, a.size()) {
      d[i] = unsigned_sat ? saturate<Lane>(a[i]) : Lane(saturate<SLane>(a[i]));
      d[a.size() + i] = unsigned_sat ? saturate<Lane>(b[i]) : Lane(saturate<SLane>(b[i]));
    }
    out.rddata = std::bit_cast<W64>(d);
    out.rdflags = 0;
  }
  return out;
}

#define sizes(b, w, d, q) ((b << 0) | (w << 1) | (d << 2) | (q << 3))

#define make_vec_implmap(name, fn, sizemask, ...)                                                                      \
  UopImpl implmap_##name[4] = {bit(sizemask, 0) ? &fn<OP_##name, W8, __VA_ARGS__> : &x86_op_nop<OP_##name, 0>,         \
                               bit(sizemask, 1) ? &fn<OP_##name, W16, __VA_ARGS__> : &x86_op_nop<OP_##name, 1>,        \
                               bit(sizemask, 2) ? &fn<OP_##name, W32, __VA_ARGS__> : &x86_op_nop<OP_##name, 2>,        \
                               bit(sizemask, 3) ? &fn<OP_##name, W64, __VA_ARGS__> : &x86_op_nop<OP_##name, 3>}

#define signed_lane(a) std::make_signed_t<decltype(a)>

make_vec_implmap(vadd, vecop, sizes(1, 1, 1, 1), [](auto a, auto b) { return decltype(a)(a + b); });
make_vec_implmap(vsub, vecop, sizes(1, 1, 1, 1), [](auto a, auto b) { return decltype(a)(a - b); });
make_vec_implmap(vadd_us, vecop, sizes(1, 1, 0, 0), [](auto a, auto b) { return saturate<decltype(a)>(W64s(a) + b); });
make_vec_implmap(vsub_us, vecop, sizes(1, 1, 0, 0), [](auto a, auto b) { return saturate<decltype(a)>(W64s(a) - b); });
make_vec_implmap(vadd_ss, vecop, sizes(1, 1, 0, 0),
                 [](auto a, auto b) { return saturate<signed_lane(a)>(W64s(signed_lane(a)(a)) + signed_lane(b)(b)); });
make_vec_implmap(vsub_ss, vecop, sizes(1, 1, 0, 0),
                 [](auto a, auto b) { return saturate<signed_lane(a)>(W64s(signed_lane(a)(a)) - signed_lane(b)(b)); });

make_vec_implmap(vshl, vecshift, sizes(0, 1, 1, 1), VecShift::Left);
make_vec_implmap(vshr, vecshift, sizes(0, 1, 1, 1), VecShift::Right);
// btv dealt with later
make_vec_implmap(vsar, vecshift, sizes(0, 1, 1, 0), VecShift::RightArith);

make_vec_implmap(vavg, vecop, sizes(1, 1, 0, 0), [](auto a, auto b) { return decltype(a)((W64(a) + b + 1) >> 1); });
// cmpv dealt with later
make_vec_implmap(vmin, vecop, sizes(1, 0, 0, 0), [](auto a, auto b) { return std::min(a, b); });
make_vec_implmap(vmax, vecop, sizes(1, 0, 0, 0), [](auto a, auto b) { return std::max(a, b); });
make_vec_implmap(vmin_s, vecop, sizes(0, 1, 0, 0),
                 [](auto a, auto b) { return std::min(signed_lane(a)(a), signed_lane(b)(b)); });
make_vec_implmap(vmax_s, vecop, sizes(0, 1, 0, 0),
                 [](auto a, auto b) { return std::max(signed_lane(a)(a), signed_lane(b)(b)); });

make_vec_implmap(vmull, vecop, sizes(0, 1, 0, 0),
                 [](auto a, auto b) { return decltype(a)(W64s(signed_lane(a)(a)) * signed_lane(b)(b)); });
// The shift count guard only placates -Wshift-count-overflow for the W64 lane
// instantiation, which is dispatched to x86_op_nop and never called.
make_vec_implmap(vmulh, vecop, sizes(0, 1, 0, 0), [](auto a, auto b) {
  return decltype(a)((W64s(signed_lane(a)(a)) * signed_lane(b)(b)) >> (sizeof(a) < 8 ? sizeof(a) * 8 : 0));
});
make_vec_implmap(vmulhu, vecop, sizes(0, 1, 0, 0),
                 [](auto a, auto b) { return decltype(a)((W64(a) * b) >> (sizeof(a) < 8 ? sizeof(a) * 8 : 0)); });

// pmaddwd: pairwise dot product of signed words into dwords. The accumulation
// is done in unsigned arithmetic: 0x8000*0x8000 + 0x8000*0x8000 wraps on x86.
UopResult uop_impl_vmaddp_w(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  const auto a = std::bit_cast<LaneArray<W16s>>(ra);
  const auto b = std::bit_cast<LaneArray<W16s>>(rb);
  LaneArray<W32> d;
  foreach (i, d.size())
    d[i] = W32(W32s(a[2 * i]) * b[2 * i]) + W32(W32s(a[2 * i + 1]) * b[2 * i + 1]);
  out.rddata = std::bit_cast<W64>(d);
  out.rdflags = 0;
  return out;
}

UopImpl implmap_vmaddp[4] = {&x86_op_nop<OP_vmaddp, 0>, &uop_impl_vmaddp_w, &x86_op_nop<OP_vmaddp, 2>,
                             &x86_op_nop<OP_vmaddp, 3>};

// psadbw: sum of absolute byte differences
UopResult uop_impl_vsad_w(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  const auto a = std::bit_cast<LaneArray<W8>>(ra);
  const auto b = std::bit_cast<LaneArray<W8>>(rb);
  W64 sum = 0;
  foreach (i, a.size())
    sum += (a[i] > b[i]) ? (a[i] - b[i]) : (b[i] - a[i]);
  out.rddata = sum;
  out.rdflags = 0;
  return out;
}

UopImpl implmap_vsad[4] = {&x86_op_nop<OP_vsad, 0>, &uop_impl_vsad_w, &x86_op_nop<OP_vsad, 2>, &x86_op_nop<OP_vsad, 3>};

make_vec_implmap(vpack_us, vecpack, sizes(1, 0, 0, 0), true);
make_vec_implmap(vpack_ss, vecpack, sizes(1, 1, 0, 0), false);

#undef signed_lane

//
// btv (bit test vector):
//
// Hardware implementation:
//
// static const W64 masks[4] = {
//   0x0101010101010101ULL,
//   0x0001000100010001ULL,
//   0x0000000100000001ULL,
//   0x0000000000000001ULL
// };
//
// int sizebits = (1 << sizeshift) * 8;
// W64 mask = masks[sizeshift] << rb;
// ra &= mask;
//
template<int sizeshift>
UopResult uop_impl_vbt(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int sizebits = (1 << sizeshift) * 8;

  rb = lowbits(rb, 3 + sizeshift);

  W64 rd = 0;

  for (int i = (1 << (3 - sizeshift)) - 1; i >= 0; i--) {
    bool b = bit(ra, (i * sizebits) + rb);
    rd = (rd << 1) | b;
  }

  out.rddata = rd;
  out.rdflags = x86_genflags<W64>(rd);
  return out;
}

UopImpl implmap_vbt[4] = {&uop_impl_vbt<0>, &uop_impl_vbt<1>, &uop_impl_vbt<2>, &uop_impl_vbt<3>};

//
// cmpv (vector compare)
// uop.cond contains the condition test
//

template<typename T>
W16 compare_and_gen_flags(T ra, T rb) {
  using U = std::make_unsigned_t<T>;

  const U a = static_cast<U>(ra);
  const U b = static_cast<U>(rb);
  const U result = U(a - b);
  const U sign_bit = U(1) << (std::numeric_limits<U>::digits - 1);

  const bool cf = a < b;
  const bool of = ((a ^ b) & (a ^ result) & sign_bit) != 0;

  return x86_genflags<T>(static_cast<T>(result)) | (cf ? FLAG_CF : 0) | (of ? FLAG_OF : 0);
}

template<int sizeshift, int cond>
UopResult uop_impl_vcmp(const UopInputs& inputs) {
  [[maybe_unused]] auto [ra, rb, rc, raflags, rbflags, rcflags, riptaken, ripseq] = inputs;
  UopResult out;
  int sizebits = (1 << sizeshift) * 8;

  W64 rd = 0;

  for (int i = (1 << (3 - sizeshift)) - 1; i >= 0; i--) {
    W64 a = bits(ra, i * sizebits, sizebits);
    W64 b = bits(rb, i * sizebits, sizebits);
    W16 flags = 0;
    switch (sizeshift) {
    case 0:
      flags = compare_and_gen_flags<byte>(a, b);
      break;
    case 1:
      flags = compare_and_gen_flags<W16>(a, b);
      break;
    case 2:
      flags = compare_and_gen_flags<W32>(a, b);
      break;
    case 3:
      flags = compare_and_gen_flags<W64>(a, b);
      break;
    }

    bool z = evaluate_cond<cond>(flags, flags);

    rd <<= sizebits;
    rd |= (z) ? bitmask(sizebits) : 0;
  }

  out.rddata = rd;
  out.rdflags = x86_genflags<W64>(rd);
  return out;
}

#define makecond(c) {&uop_impl_vcmp<0, c>, &uop_impl_vcmp<1, c>, &uop_impl_vcmp<2, c>, &uop_impl_vcmp<3, c>}

UopImpl implmap_vcmp[16][4] = {makecond(0),  makecond(1),  makecond(2),  makecond(3), makecond(4),  makecond(5),
                               makecond(6),  makecond(7),  makecond(8),  makecond(9), makecond(10), makecond(11),
                               makecond(12), makecond(13), makecond(14), makecond(15)};

#undef makecond
#undef sizes

UopImpl get_synthcode_for_uop(int op, int size, bool setflags, int cond, int extshift, bool except, bool internal) {
  UopImpl func;

  switch (op) {
  case OP_nop:
    func = &uop_impl_nop;
    break;
  case OP_mov:
    func = implmap_mov[size];
    break;
  case OP_and:
    func = implmap_and[size][setflags];
    break;
  case OP_or:
    func = implmap_or[size][setflags];
    break;
  case OP_xor:
    func = implmap_xor[size][setflags];
    break;
  case OP_andnot:
    func = implmap_andnot[size][setflags];
    break;
  case OP_ornot:
    func = implmap_ornot[size][setflags];
    break;
  case OP_nand:
    func = implmap_nand[size][setflags];
    break;
  case OP_nor:
    func = implmap_nor[size][setflags];
    break;
  case OP_eqv:
    func = implmap_eqv[size][setflags];
    break;
  case OP_add:
    func = implmap_add[size][setflags];
    break;
  case OP_sub:
    func = implmap_sub[size][setflags];
    break;
  case OP_adda:
    func = implmap_adda[size][extshift][setflags];
    break;
  case OP_suba:
    func = implmap_suba[size][extshift][setflags];
    break;
  case OP_addm:
    func = implmap_addm[size][setflags];
    break;
  case OP_subm:
    func = implmap_subm[size][setflags];
    break;
  case OP_sel:
    func = implmap_sel[cond][size];
    break;
  case OP_sel_cmp:
    func = implmap_sel_cmp[cond][size][extshift];
    break;
  case OP_set:
    func = implmap_set[cond][size];
    break;
  case OP_set_and:
    func = implmap_set_and[cond][size][extshift];
    break;
  case OP_set_sub:
    func = implmap_set_sub[cond][size][extshift];
    break;
  case OP_br:
    func = implmap_br[cond][except];
    break;
  case OP_br_sub:
    func = implmap_br_sub[cond][size][except];
    break;
  case OP_br_and:
    func = implmap_br_and[cond][size][except];
    break;
  case OP_jmp:
    func = (except ? &uop_impl_jmp_ex : &uop_impl_jmp);
    break;
  case OP_bru:
    func = &uop_impl_bru;
    break;
  case OP_brp:
    func = &uop_impl_brp;
    break;
  case OP_chk:
    func = implmap_chk[cond][0];
    break;
  case OP_chk_sub:
    func = implmap_chk_sub[cond][size];
    break;
  case OP_chk_and:
    func = implmap_chk_and[cond][size];
    break;

    //
    // Loads and stores are handled specially in the core model:
    //
  case OP_ld:
  case OP_ldx:
  case OP_ld_pre:
  case OP_ld_a16:
  case OP_st:
  case OP_st_a16:
  case OP_mf:
    func = &uop_impl_nop;
    break;

  case OP_bt:
    func = implmap_bt[size][setflags];
    break;
  case OP_bts:
    func = implmap_bts[size][setflags];
    break;
  case OP_btr:
    func = implmap_btr[size][setflags];
    break;
  case OP_btc:
    func = implmap_btc[size][setflags];
    break;

  case OP_rotl:
    func = implmap_rotl[size][setflags];
    break;
  case OP_rotr:
    func = implmap_rotr[size][setflags];
    break;
  case OP_rotcl:
    func = implmap_rotcl[size][setflags];
    break;
  case OP_rotcr:
    func = implmap_rotcr[size][setflags];
    break;
  case OP_shl:
    func = implmap_shl[size][setflags];
    break;
  case OP_shr:
    func = implmap_shr[size][setflags];
    break;
  case OP_sar:
    func = implmap_sar[size][setflags];
    break;
  case OP_mask:
    func = implmap_mask[size][cond];
    break;

  case OP_shls:
    func = implmap_shls[size][setflags];
    break;
  case OP_shrs:
    func = implmap_shrs[size][setflags];
    break;
  case OP_sars:
    func = implmap_sars[size][setflags];
    break;
  case OP_maskb:
    func = implmap_maskb[size][cond];
    break;

  case OP_bswap:
    func = implmap_bswap[size][0];
    break;
  case OP_collcc:
    func = &uop_impl_collcc;
    break;
  case OP_movccr:
    func = &uop_impl_movccr;
    break;
  case OP_movrcc:
    func = &uop_impl_movrcc;
    break;
  case OP_andcc:
    func = &uop_impl_andcc;
    break;
  case OP_orcc:
    func = &uop_impl_orcc;
    break;
  case OP_ornotcc:
    func = &uop_impl_ornotcc;
    break;
  case OP_xorcc:
    func = &uop_impl_xorcc;
    break;

  case OP_mull:
    func = implmap_mull[size][setflags];
    break;
  case OP_mulh:
    func = implmap_mulh[size][setflags];
    break;
  case OP_mulhu:
    func = implmap_mulhu[size][setflags];
    break;
  case OP_mulhl:
    func = implmap_mulhl[size];
    break;

  case OP_ctz:
    func = implmap_ctz[size][setflags];
    break;
  case OP_clz:
    func = implmap_clz[size][setflags];
    break;
    // case OP_ctpop:
  case OP_permb:
    func = &uop_impl_permb;
    break;

  case OP_div:
    func = implmap_div[size];
    break;
  case OP_rem:
    func = implmap_rem[size];
    break;
  case OP_divs:
    func = implmap_divs[size];
    break;
  case OP_rems:
    func = implmap_rems[size];
    break;

  case OP_min:
    func = implmap_min[size];
    break;
  case OP_max:
    func = implmap_max[size];
    break;
  case OP_min_s:
    func = implmap_min_s[size];
    break;
  case OP_max_s:
    func = implmap_max_s[size];
    break;

  case OP_fadd:
    func = implmap_fadd[size];
    break;
  case OP_fsub:
    func = implmap_fsub[size];
    break;
  case OP_fmul:
    func = implmap_fmul[size];
    break;
  case OP_fmadd:
    func = implmap_fmadd[size];
    break;
  case OP_fmsub:
    func = implmap_fmsub[size];
    break;
  case OP_fmsubr:
    func = implmap_fmsubr[size];
    break;
  case OP_fdiv:
    func = implmap_fdiv[size];
    break;
  case OP_fsqrt:
    func = implmap_fsqrt[size];
    break;
  case OP_frcp:
    func = implmap_frcp[size];
    break;
  case OP_frsqrt:
    func = implmap_frsqrt[size];
    break;
  case OP_fmin:
    func = implmap_fmin[size];
    break;
  case OP_fmax:
    func = implmap_fmax[size];
    break;
  case OP_fcmp:
    func = implmap_fcmp[cond][size];
    break;
  case OP_fcmpcc:
    func = implmap_fcmpcc[cond][size];
    break;

  case OP_fcvt_i2s_ins:
    func = &uop_impl_fcvt_i2s_ins;
    break;

  case OP_fcvt_i2s_p:
    func = &uop_impl_fcvt_i2s_p;
    break;
  case OP_fcvt_i2d_lo:
    func = &uop_impl_fcvt_i2d_lo;
    break;
  case OP_fcvt_i2d_hi:
    func = &uop_impl_fcvt_i2d_hi;
    break;

  case OP_fcvt_q2s_ins:
    func = &uop_impl_fcvt_q2s_ins;
    break;
  case OP_fcvt_q2d:
    func = &uop_impl_fcvt_q2d;
    break;

  case OP_fcvt_s2i:
    func = implmap_fcvt_s2i[size & 1];
    break;
  case OP_fcvt_s2q:
    func = implmap_fcvt_s2q[size & 1];
    break;
  case OP_fcvt_s2i_p:
    func = implmap_fcvt_s2i_p[size & 1];
    break;
  case OP_fcvt_d2i:
    func = implmap_fcvt_d2i[size & 1];
    break;
  case OP_fcvt_d2q:
    func = implmap_fcvt_d2q[size & 1];
    break;
  case OP_fcvt_d2i_p:
    func = implmap_fcvt_d2i_p[size & 1];
    break;
  case OP_fcvt_d2s_ins:
    func = &uop_impl_fcvt_d2s_ins;
    break;
  case OP_fcvt_d2s_p:
    func = implmap_fcvt_d2s_p[0];
    break;

  case OP_fcvt_s2d_lo:
    func = &uop_impl_fcvt_s2d_lo;
    break;
  case OP_fcvt_s2d_hi:
    func = &uop_impl_fcvt_s2d_hi;
    break;

  case OP_vadd:
    func = implmap_vadd[size];
    break;
  case OP_vsub:
    func = implmap_vsub[size];
    break;
  case OP_vadd_us:
    func = implmap_vadd_us[size];
    break;
  case OP_vsub_us:
    func = implmap_vsub_us[size];
    break;
  case OP_vadd_ss:
    func = implmap_vadd_ss[size];
    break;
  case OP_vsub_ss:
    func = implmap_vsub_ss[size];
    break;
  case OP_vshl:
    func = implmap_vshl[size];
    break;
  case OP_vshr:
    func = implmap_vshr[size];
    break;
  case OP_vbt:
    func = implmap_vbt[size];
    break;
  case OP_vsar:
    func = implmap_vsar[size];
    break;
  case OP_vavg:
    func = implmap_vavg[size];
    break;
  case OP_vcmp:
    func = implmap_vcmp[cond][size];
    break;
  case OP_vmin:
    func = implmap_vmin[size];
    break;
  case OP_vmax:
    func = implmap_vmax[size];
    break;
  case OP_vmin_s:
    func = implmap_vmin_s[size];
    break;
  case OP_vmax_s:
    func = implmap_vmax_s[size];
    break;
  case OP_vmull:
    func = implmap_vmull[size];
    break;
  case OP_vmulh:
    func = implmap_vmulh[size];
    break;
  case OP_vmulhu:
    func = implmap_vmulhu[size];
    break;
  case OP_vmaddp:
    func = implmap_vmaddp[size];
    break;
  case OP_vsad:
    func = implmap_vsad[size];
    break;
  case OP_vpack_us:
    func = implmap_vpack_us[size];
    break;
  case OP_vpack_ss:
    func = implmap_vpack_ss[size];
    break;
  default:
    logging::println("Unknown uop opcode {} ({})", op, nameof(op));
    logging::flush();
    assert(false);
  }
  return func;
}

void synth_uops_for_bb(BasicBlock& bb) {
  bb.synthops = new UopImpl[bb.count];
  foreach (i, bb.count) {
    const TransOp& transop = bb.transops[i];
    UopImpl func = get_synthcode_for_uop(transop.opcode, transop.size, transop.setflags, transop.cond, transop.extshift,
                                         0, transop.internal);
    bb.synthops[i] = func;
  }
}

UopImpl get_synthcode_for_cond_branch(int opcode, int cond, int size, bool except) {
  UopImpl func;

  switch (opcode) {
  case OP_br_sub:
    func = implmap_br_sub[cond][size][except];
    break;
  case OP_br_and:
    func = implmap_br_and[cond][size][except];
    break;
  case OP_br:
    func = implmap_br[cond][except];
    break;
  default:
    assert(false);
  }

  return func;
}

void init_uops() {}

void shutdown_uops() {}

} // namespace x86sim
