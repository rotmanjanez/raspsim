//
// Super Standard Template Library
//
// Copyright 1997-2008 Matt T. Yourst <yourst@yourst.com>
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
//

#include "globals.h"
#include "superstl.h"

namespace x86sim {

// For debugging of messages before crashes:
bool force_synchronous_streams = false;

namespace superstl {
//
// x86 compatible non-excepting divide and remainder:
//
// This code is heavily modified from the Bochs version:
//
struct W128 {
  W64 lo;
  W64 hi;
};

struct W128s {
  W64 lo;
  W64s hi;
};

void long_neg(W128s& n) {
  W64 t = n.lo;
  n.lo = -n.lo;
  if (t - 1 > t)
    --n.hi;
  n.hi = ~n.hi;
}

void long_shl(W128& a) {
  W64 c;
  c = a.lo >> 63;
  a.lo <<= 1;
  a.hi <<= 1;
  a.hi |= c;
}

void long_shr(W128& a) {
  W64 c;
  c = a.hi << 63;
  a.hi >>= 1;
  a.lo >>= 1;
  a.lo |= c;
}

bool long_sub(W128& a, W128& b) {
  W64 t = a.lo;
  a.lo -= b.lo;
  int c = (a.lo > t);
  t = a.hi;
  a.hi -= b.hi + c;
  return (a.hi > t);
}

bool long_le(W128& a, W128& b) {
  if (a.hi == b.hi) {
    return (a.lo <= b.lo);
  } else {
    return (a.hi <= b.hi);
  }
}

void long_div(W128& quotient, W64& remainder, const W128& dividend, W64 divisor) {
  W128 d, acc, q, temp;
  int n, c;

  d.lo = divisor;
  d.hi = 0;
  acc.lo = dividend.lo;
  acc.hi = dividend.hi;
  q.lo = 0;
  q.hi = 0;
  n = 0;

  while (long_le(d, acc) && (n < 128)) {
    long_shl(d);
    n++;
  }

  while (n > 0) {
    long_shr(d);
    long_shl(q);
    temp.lo = acc.lo;
    temp.hi = acc.hi;
    c = long_sub(acc, d);
    if (c) {
      acc.lo = temp.lo;
      acc.hi = temp.hi;
    } else {
      q.lo++;
    }
    n--;
  }

  remainder = acc.lo;
  quotient.lo = q.lo;
  quotient.hi = q.hi;
}

void long_idiv(W128s& quotient, W64s& remainder, W128s& dividend, W64s divisor) {
  W128s temp = dividend;

  bool dividend_was_negative = (temp.hi < 0);
  bool divisor_was_negative = (divisor < 0);
  if (dividend_was_negative)
    long_neg(temp);
  if (divisor_was_negative)
    divisor = -divisor;

  long_div((W128&)quotient, (W64&)remainder, (W128&)temp, divisor);

  if (dividend_was_negative ^ divisor_was_negative)
    long_neg(quotient);

  // NOTE: Bochs original code was: if (divisor_was_negative) remainder = -remainder;
  // This is actually a bug: Intel and AMD manuals say sign of remainder is equal to sign of dividend, not divisor.
  if (dividend_was_negative)
    remainder = -remainder;
}

template<>
bool div_rem(W64& quotientlo, W64& remainder, W64 dividend_hi, W64 dividend_lo, W64 divisor) {
  W128 dividend;
  W128 quotient;

  dividend.lo = dividend_lo;
  dividend.hi = dividend_hi;

  if unlikely (!divisor)
    goto out;

  long_div(quotient, remainder, dividend, divisor);
  quotientlo = quotient.lo;

  if unlikely (quotient.hi != 0)
    goto out;

  return true;

out:
  quotientlo = 0;
  remainder = 0;
  return false;
}

template<>
bool div_rem_s(W64& quotientlo, W64& remainder, W64 dividend_hi, W64 dividend_lo, W64 divisor) {
  W64s op2_64, remainder_64, quotient_64l;
  W128s dividend;
  W128s quotient;

  dividend.lo = dividend_lo;
  dividend.hi = dividend_hi;

  if unlikely (!divisor)
    goto out;

  if unlikely ((divisor == -1) && (dividend.hi == 0x8000000000000000ULL) && (!dividend.lo))
    goto out;

  long_idiv(quotient, (W64s&)remainder, dividend, divisor);
  quotientlo = quotient.lo;

  if unlikely ((!(quotient.lo & 0x8000000000000000ULL) && quotient.hi != 0) ||
               (quotient.lo & 0x8000000000000000ULL) && quotient.hi != 0xffffffffffffffffULL)
    goto out;

  return true;

out:
  quotientlo = 0;
  remainder = 0;
  return false;
}

template<typename T>
bool div_rem(T& quotientlo, T& remainder, T dividend_hi, T dividend_lo, T divisor) {
  static const int B = (sizeof(T) * 8);

  W64 dividend = (W64(dividend_hi) << B) + dividend_lo;
  W64 quotient;

  if unlikely (divisor == 0)
    goto out;

  quotient = dividend / divisor;
  remainder = T(dividend % divisor);
  quotientlo = T(quotient);

  if unlikely (quotient != quotientlo)
    goto out;

  return true;

out:
  quotientlo = 0;
  remainder = 0;
  return false;
}

template<typename T>
bool div_rem_s(T& quotientlo, T& remainder, T dividend_hi, T dividend_lo, T udivisor) {
  static const int B = (sizeof(T) * 8);

  W64s dividend = (W64(signext64(dividend_hi, B)) << B) + dividend_lo;
  W64 quotient;
  W64s divisor = signext64(W64(udivisor), B);

  if unlikely (divisor == 0)
    goto out;

  // check MIN_INT divided by -1 case
  // e.g. MIN_INT for 32/16/8 = 0x8000000000000000ULL, 0x80000000ULL, 0x8000ULL
  if unlikely ((dividend == signext64((1ULL << ((B * 2) - 1)), (B * 2))) && (divisor == bitmask(B)))
    goto out;

  quotient = dividend / divisor;
  remainder = T(dividend % divisor);
  quotientlo = T(quotient);

  if unlikely (quotient != signext64(W64(quotientlo), B))
    goto out;

  return true;

out:
  quotientlo = 0;
  remainder = 0;
  return false;
}

#define decl_div_rem(T) template bool div_rem(T& quotient, T& remainder, T dividend_hi, T dividend_lo, T divisor);
#define decl_div_rem_s(T) template bool div_rem_s(T& quotient, T& remainder, T dividend_hi, T dividend_lo, T divisor);

decl_div_rem(W8);
decl_div_rem(W16);
decl_div_rem(W32);

decl_div_rem_s(W8);
decl_div_rem_s(W16);
decl_div_rem_s(W32);

}; // namespace superstl

using namespace superstl;

const unsigned char popcountlut8bit[] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2,
    3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3,
    3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
    6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4,
    3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4,
    5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
    6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

const byte lsbindexlut8bit[256] = {
    0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2,
    0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0,
    1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
    2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
    1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};

const W32 CRC32::crctable[256] = {
    0x00000000L, 0x77073096L, 0xee0e612cL, 0x990951baL, 0x076dc419L, 0x706af48fL, 0xe963a535L, 0x9e6495a3L, 0x0edb8832L,
    0x79dcb8a4L, 0xe0d5e91eL, 0x97d2d988L, 0x09b64c2bL, 0x7eb17cbdL, 0xe7b82d07L, 0x90bf1d91L, 0x1db71064L, 0x6ab020f2L,
    0xf3b97148L, 0x84be41deL, 0x1adad47dL, 0x6ddde4ebL, 0xf4d4b551L, 0x83d385c7L, 0x136c9856L, 0x646ba8c0L, 0xfd62f97aL,
    0x8a65c9ecL, 0x14015c4fL, 0x63066cd9L, 0xfa0f3d63L, 0x8d080df5L, 0x3b6e20c8L, 0x4c69105eL, 0xd56041e4L, 0xa2677172L,
    0x3c03e4d1L, 0x4b04d447L, 0xd20d85fdL, 0xa50ab56bL, 0x35b5a8faL, 0x42b2986cL, 0xdbbbc9d6L, 0xacbcf940L, 0x32d86ce3L,
    0x45df5c75L, 0xdcd60dcfL, 0xabd13d59L, 0x26d930acL, 0x51de003aL, 0xc8d75180L, 0xbfd06116L, 0x21b4f4b5L, 0x56b3c423L,
    0xcfba9599L, 0xb8bda50fL, 0x2802b89eL, 0x5f058808L, 0xc60cd9b2L, 0xb10be924L, 0x2f6f7c87L, 0x58684c11L, 0xc1611dabL,
    0xb6662d3dL, 0x76dc4190L, 0x01db7106L, 0x98d220bcL, 0xefd5102aL, 0x71b18589L, 0x06b6b51fL, 0x9fbfe4a5L, 0xe8b8d433L,
    0x7807c9a2L, 0x0f00f934L, 0x9609a88eL, 0xe10e9818L, 0x7f6a0dbbL, 0x086d3d2dL, 0x91646c97L, 0xe6635c01L, 0x6b6b51f4L,
    0x1c6c6162L, 0x856530d8L, 0xf262004eL, 0x6c0695edL, 0x1b01a57bL, 0x8208f4c1L, 0xf50fc457L, 0x65b0d9c6L, 0x12b7e950L,
    0x8bbeb8eaL, 0xfcb9887cL, 0x62dd1ddfL, 0x15da2d49L, 0x8cd37cf3L, 0xfbd44c65L, 0x4db26158L, 0x3ab551ceL, 0xa3bc0074L,
    0xd4bb30e2L, 0x4adfa541L, 0x3dd895d7L, 0xa4d1c46dL, 0xd3d6f4fbL, 0x4369e96aL, 0x346ed9fcL, 0xad678846L, 0xda60b8d0L,
    0x44042d73L, 0x33031de5L, 0xaa0a4c5fL, 0xdd0d7cc9L, 0x5005713cL, 0x270241aaL, 0xbe0b1010L, 0xc90c2086L, 0x5768b525L,
    0x206f85b3L, 0xb966d409L, 0xce61e49fL, 0x5edef90eL, 0x29d9c998L, 0xb0d09822L, 0xc7d7a8b4L, 0x59b33d17L, 0x2eb40d81L,
    0xb7bd5c3bL, 0xc0ba6cadL, 0xedb88320L, 0x9abfb3b6L, 0x03b6e20cL, 0x74b1d29aL, 0xead54739L, 0x9dd277afL, 0x04db2615L,
    0x73dc1683L, 0xe3630b12L, 0x94643b84L, 0x0d6d6a3eL, 0x7a6a5aa8L, 0xe40ecf0bL, 0x9309ff9dL, 0x0a00ae27L, 0x7d079eb1L,
    0xf00f9344L, 0x8708a3d2L, 0x1e01f268L, 0x6906c2feL, 0xf762575dL, 0x806567cbL, 0x196c3671L, 0x6e6b06e7L, 0xfed41b76L,
    0x89d32be0L, 0x10da7a5aL, 0x67dd4accL, 0xf9b9df6fL, 0x8ebeeff9L, 0x17b7be43L, 0x60b08ed5L, 0xd6d6a3e8L, 0xa1d1937eL,
    0x38d8c2c4L, 0x4fdff252L, 0xd1bb67f1L, 0xa6bc5767L, 0x3fb506ddL, 0x48b2364bL, 0xd80d2bdaL, 0xaf0a1b4cL, 0x36034af6L,
    0x41047a60L, 0xdf60efc3L, 0xa867df55L, 0x316e8eefL, 0x4669be79L, 0xcb61b38cL, 0xbc66831aL, 0x256fd2a0L, 0x5268e236L,
    0xcc0c7795L, 0xbb0b4703L, 0x220216b9L, 0x5505262fL, 0xc5ba3bbeL, 0xb2bd0b28L, 0x2bb45a92L, 0x5cb36a04L, 0xc2d7ffa7L,
    0xb5d0cf31L, 0x2cd99e8bL, 0x5bdeae1dL, 0x9b64c2b0L, 0xec63f226L, 0x756aa39cL, 0x026d930aL, 0x9c0906a9L, 0xeb0e363fL,
    0x72076785L, 0x05005713L, 0x95bf4a82L, 0xe2b87a14L, 0x7bb12baeL, 0x0cb61b38L, 0x92d28e9bL, 0xe5d5be0dL, 0x7cdcefb7L,
    0x0bdbdf21L, 0x86d3d2d4L, 0xf1d4e242L, 0x68ddb3f8L, 0x1fda836eL, 0x81be16cdL, 0xf6b9265bL, 0x6fb077e1L, 0x18b74777L,
    0x88085ae6L, 0xff0f6a70L, 0x66063bcaL, 0x11010b5cL, 0x8f659effL, 0xf862ae69L, 0x616bffd3L, 0x166ccf45L, 0xa00ae278L,
    0xd70dd2eeL, 0x4e048354L, 0x3903b3c2L, 0xa7672661L, 0xd06016f7L, 0x4969474dL, 0x3e6e77dbL, 0xaed16a4aL, 0xd9d65adcL,
    0x40df0b66L, 0x37d83bf0L, 0xa9bcae53L, 0xdebb9ec5L, 0x47b2cf7fL, 0x30b5ffe9L, 0xbdbdf21cL, 0xcabac28aL, 0x53b39330L,
    0x24b4a3a6L, 0xbad03605L, 0xcdd70693L, 0x54de5729L, 0x23d967bfL, 0xb3667a2eL, 0xc4614ab8L, 0x5d681b02L, 0x2a6f2b94L,
    0xb40bbe37L, 0xc30c8ea1L, 0x5a05df1bL, 0x2d02ef8dL};

const W64 expand_8bit_to_64bit_lut[256] alignto(8) = {
    0x0000000000000000ULL, 0x00000000000000ffULL, 0x000000000000ff00ULL, 0x000000000000ffffULL, 0x0000000000ff0000ULL,
    0x0000000000ff00ffULL, 0x0000000000ffff00ULL, 0x0000000000ffffffULL, 0x00000000ff000000ULL, 0x00000000ff0000ffULL,
    0x00000000ff00ff00ULL, 0x00000000ff00ffffULL, 0x00000000ffff0000ULL, 0x00000000ffff00ffULL, 0x00000000ffffff00ULL,
    0x00000000ffffffffULL, 0x000000ff00000000ULL, 0x000000ff000000ffULL, 0x000000ff0000ff00ULL, 0x000000ff0000ffffULL,
    0x000000ff00ff0000ULL, 0x000000ff00ff00ffULL, 0x000000ff00ffff00ULL, 0x000000ff00ffffffULL, 0x000000ffff000000ULL,
    0x000000ffff0000ffULL, 0x000000ffff00ff00ULL, 0x000000ffff00ffffULL, 0x000000ffffff0000ULL, 0x000000ffffff00ffULL,
    0x000000ffffffff00ULL, 0x000000ffffffffffULL, 0x0000ff0000000000ULL, 0x0000ff00000000ffULL, 0x0000ff000000ff00ULL,
    0x0000ff000000ffffULL, 0x0000ff0000ff0000ULL, 0x0000ff0000ff00ffULL, 0x0000ff0000ffff00ULL, 0x0000ff0000ffffffULL,
    0x0000ff00ff000000ULL, 0x0000ff00ff0000ffULL, 0x0000ff00ff00ff00ULL, 0x0000ff00ff00ffffULL, 0x0000ff00ffff0000ULL,
    0x0000ff00ffff00ffULL, 0x0000ff00ffffff00ULL, 0x0000ff00ffffffffULL, 0x0000ffff00000000ULL, 0x0000ffff000000ffULL,
    0x0000ffff0000ff00ULL, 0x0000ffff0000ffffULL, 0x0000ffff00ff0000ULL, 0x0000ffff00ff00ffULL, 0x0000ffff00ffff00ULL,
    0x0000ffff00ffffffULL, 0x0000ffffff000000ULL, 0x0000ffffff0000ffULL, 0x0000ffffff00ff00ULL, 0x0000ffffff00ffffULL,
    0x0000ffffffff0000ULL, 0x0000ffffffff00ffULL, 0x0000ffffffffff00ULL, 0x0000ffffffffffffULL, 0x00ff000000000000ULL,
    0x00ff0000000000ffULL, 0x00ff00000000ff00ULL, 0x00ff00000000ffffULL, 0x00ff000000ff0000ULL, 0x00ff000000ff00ffULL,
    0x00ff000000ffff00ULL, 0x00ff000000ffffffULL, 0x00ff0000ff000000ULL, 0x00ff0000ff0000ffULL, 0x00ff0000ff00ff00ULL,
    0x00ff0000ff00ffffULL, 0x00ff0000ffff0000ULL, 0x00ff0000ffff00ffULL, 0x00ff0000ffffff00ULL, 0x00ff0000ffffffffULL,
    0x00ff00ff00000000ULL, 0x00ff00ff000000ffULL, 0x00ff00ff0000ff00ULL, 0x00ff00ff0000ffffULL, 0x00ff00ff00ff0000ULL,
    0x00ff00ff00ff00ffULL, 0x00ff00ff00ffff00ULL, 0x00ff00ff00ffffffULL, 0x00ff00ffff000000ULL, 0x00ff00ffff0000ffULL,
    0x00ff00ffff00ff00ULL, 0x00ff00ffff00ffffULL, 0x00ff00ffffff0000ULL, 0x00ff00ffffff00ffULL, 0x00ff00ffffffff00ULL,
    0x00ff00ffffffffffULL, 0x00ffff0000000000ULL, 0x00ffff00000000ffULL, 0x00ffff000000ff00ULL, 0x00ffff000000ffffULL,
    0x00ffff0000ff0000ULL, 0x00ffff0000ff00ffULL, 0x00ffff0000ffff00ULL, 0x00ffff0000ffffffULL, 0x00ffff00ff000000ULL,
    0x00ffff00ff0000ffULL, 0x00ffff00ff00ff00ULL, 0x00ffff00ff00ffffULL, 0x00ffff00ffff0000ULL, 0x00ffff00ffff00ffULL,
    0x00ffff00ffffff00ULL, 0x00ffff00ffffffffULL, 0x00ffffff00000000ULL, 0x00ffffff000000ffULL, 0x00ffffff0000ff00ULL,
    0x00ffffff0000ffffULL, 0x00ffffff00ff0000ULL, 0x00ffffff00ff00ffULL, 0x00ffffff00ffff00ULL, 0x00ffffff00ffffffULL,
    0x00ffffffff000000ULL, 0x00ffffffff0000ffULL, 0x00ffffffff00ff00ULL, 0x00ffffffff00ffffULL, 0x00ffffffffff0000ULL,
    0x00ffffffffff00ffULL, 0x00ffffffffffff00ULL, 0x00ffffffffffffffULL, 0xff00000000000000ULL, 0xff000000000000ffULL,
    0xff0000000000ff00ULL, 0xff0000000000ffffULL, 0xff00000000ff0000ULL, 0xff00000000ff00ffULL, 0xff00000000ffff00ULL,
    0xff00000000ffffffULL, 0xff000000ff000000ULL, 0xff000000ff0000ffULL, 0xff000000ff00ff00ULL, 0xff000000ff00ffffULL,
    0xff000000ffff0000ULL, 0xff000000ffff00ffULL, 0xff000000ffffff00ULL, 0xff000000ffffffffULL, 0xff0000ff00000000ULL,
    0xff0000ff000000ffULL, 0xff0000ff0000ff00ULL, 0xff0000ff0000ffffULL, 0xff0000ff00ff0000ULL, 0xff0000ff00ff00ffULL,
    0xff0000ff00ffff00ULL, 0xff0000ff00ffffffULL, 0xff0000ffff000000ULL, 0xff0000ffff0000ffULL, 0xff0000ffff00ff00ULL,
    0xff0000ffff00ffffULL, 0xff0000ffffff0000ULL, 0xff0000ffffff00ffULL, 0xff0000ffffffff00ULL, 0xff0000ffffffffffULL,
    0xff00ff0000000000ULL, 0xff00ff00000000ffULL, 0xff00ff000000ff00ULL, 0xff00ff000000ffffULL, 0xff00ff0000ff0000ULL,
    0xff00ff0000ff00ffULL, 0xff00ff0000ffff00ULL, 0xff00ff0000ffffffULL, 0xff00ff00ff000000ULL, 0xff00ff00ff0000ffULL,
    0xff00ff00ff00ff00ULL, 0xff00ff00ff00ffffULL, 0xff00ff00ffff0000ULL, 0xff00ff00ffff00ffULL, 0xff00ff00ffffff00ULL,
    0xff00ff00ffffffffULL, 0xff00ffff00000000ULL, 0xff00ffff000000ffULL, 0xff00ffff0000ff00ULL, 0xff00ffff0000ffffULL,
    0xff00ffff00ff0000ULL, 0xff00ffff00ff00ffULL, 0xff00ffff00ffff00ULL, 0xff00ffff00ffffffULL, 0xff00ffffff000000ULL,
    0xff00ffffff0000ffULL, 0xff00ffffff00ff00ULL, 0xff00ffffff00ffffULL, 0xff00ffffffff0000ULL, 0xff00ffffffff00ffULL,
    0xff00ffffffffff00ULL, 0xff00ffffffffffffULL, 0xffff000000000000ULL, 0xffff0000000000ffULL, 0xffff00000000ff00ULL,
    0xffff00000000ffffULL, 0xffff000000ff0000ULL, 0xffff000000ff00ffULL, 0xffff000000ffff00ULL, 0xffff000000ffffffULL,
    0xffff0000ff000000ULL, 0xffff0000ff0000ffULL, 0xffff0000ff00ff00ULL, 0xffff0000ff00ffffULL, 0xffff0000ffff0000ULL,
    0xffff0000ffff00ffULL, 0xffff0000ffffff00ULL, 0xffff0000ffffffffULL, 0xffff00ff00000000ULL, 0xffff00ff000000ffULL,
    0xffff00ff0000ff00ULL, 0xffff00ff0000ffffULL, 0xffff00ff00ff0000ULL, 0xffff00ff00ff00ffULL, 0xffff00ff00ffff00ULL,
    0xffff00ff00ffffffULL, 0xffff00ffff000000ULL, 0xffff00ffff0000ffULL, 0xffff00ffff00ff00ULL, 0xffff00ffff00ffffULL,
    0xffff00ffffff0000ULL, 0xffff00ffffff00ffULL, 0xffff00ffffffff00ULL, 0xffff00ffffffffffULL, 0xffffff0000000000ULL,
    0xffffff00000000ffULL, 0xffffff000000ff00ULL, 0xffffff000000ffffULL, 0xffffff0000ff0000ULL, 0xffffff0000ff00ffULL,
    0xffffff0000ffff00ULL, 0xffffff0000ffffffULL, 0xffffff00ff000000ULL, 0xffffff00ff0000ffULL, 0xffffff00ff00ff00ULL,
    0xffffff00ff00ffffULL, 0xffffff00ffff0000ULL, 0xffffff00ffff00ffULL, 0xffffff00ffffff00ULL, 0xffffff00ffffffffULL,
    0xffffffff00000000ULL, 0xffffffff000000ffULL, 0xffffffff0000ff00ULL, 0xffffffff0000ffffULL, 0xffffffff00ff0000ULL,
    0xffffffff00ff00ffULL, 0xffffffff00ffff00ULL, 0xffffffff00ffffffULL, 0xffffffffff000000ULL, 0xffffffffff0000ffULL,
    0xffffffffff00ff00ULL, 0xffffffffff00ffffULL, 0xffffffffffff0000ULL, 0xffffffffffff00ffULL, 0xffffffffffffff00ULL,
    0xffffffffffffffffULL,
};

} // namespace x86sim
