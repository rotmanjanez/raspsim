// -*- c++ -*-
//
// Cross-platform globals - refactored for C++23
//
// Originally from PTLsim, Copyright 1997-2008 Matt T. Yourst <yourst@yourst.com>
// Refactored for architecture independence
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
//

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <bit>
#include <atomic>
#include <chrono>
#include <limits>
#include <type_traits>

extern "C" {
#include <sys/ptrace.h>
}

#include <typedefs.h>

#ifdef __cplusplus

#include <cfloat>

#define __stringify_1(x) #x
#define stringify(x) __stringify_1(x)

#define alignto(x) __attribute__ ((aligned (x)))
#define insection(x) __attribute__ ((section (x)))
#define packedstruct __attribute__ ((packed))
#define noinline __attribute__((noinline))

#define unlikely(x) (__builtin_expect(!!(x), 0))
#define likely(x) (__builtin_expect(!!(x), 1))
#define isconst(x) (__builtin_constant_p(x))
#define getcaller() (__builtin_return_address(0))
#define asmlinkage extern "C"

//
// Asserts
//
#if defined __cplusplus
#  define __ASSERT_VOID_CAST static_cast<void>
#else
#  define __ASSERT_VOID_CAST (void)
#endif

asmlinkage void assert_fail(const char *__assertion, const char *__file, unsigned int __line, const char *__function) __attribute__ ((__noreturn__));

// For embedded debugging use only - portable version:
static inline void assert_fail_trap(const char *__assertion, const char *__file, unsigned int __line, const char *__function) {
  std::abort();
}

#define __CONCAT(x,y)	x ## y
#define __STRING(x)	#x

#define nan NAN
#define inf INFINITY

template <typename T> struct limits { static const T min = 0; static const T max = 0; };
#define MakeLimits(T, __min, __max) template <> struct limits<T> { static const T min = (__min); static const T max = (__max); };
MakeLimits(W8, 0, 0xff);
MakeLimits(W16, 0, 0xffff);
MakeLimits(W32, 0, 0xffffffff);
MakeLimits(W64, 0, 0xffffffffffffffffULL);
MakeLimits(W8s, 0x80, 0x7f);
MakeLimits(W16s, 0x8000, 0x7fff);
MakeLimits(W32s, 0x80000000, 0x7fffffff);
MakeLimits(W64s, 0x8000000000000000LL, 0x7fffffffffffffffLL);
// Note: On 64-bit systems, long/unsigned long are same type as W64s/W64
#undef MakeLimits

template <typename T> struct isprimitive_t { static const bool primitive = 0; };
#define MakePrimitive(T) template <> struct isprimitive_t<T> { static const bool primitive = 1; }
MakePrimitive(signed char);
MakePrimitive(unsigned char);
MakePrimitive(signed short);
MakePrimitive(unsigned short);
MakePrimitive(signed int);
MakePrimitive(unsigned int);
MakePrimitive(signed long);
MakePrimitive(unsigned long);
MakePrimitive(signed long long);
MakePrimitive(unsigned long long);
MakePrimitive(float);
MakePrimitive(double);
MakePrimitive(bool);

template<typename T> struct ispointer_t { static const bool pointer = 0; };
template <typename T> struct ispointer_t<T*> { static const bool pointer = 1; };
#define ispointer(T) (ispointer_t<T>::pointer)
#define isprimitive(T) (isprimitive_t<T>::primitive)

// Null pointer to the specified object type, for computing field offsets
#define offsetof_(T, field) ((Waddr)(&(reinterpret_cast<T*>(0)->field)) - ((Waddr)reinterpret_cast<T*>(0)))
#define baseof(T, field, ptr) ((T*)(((byte*)(ptr)) - offsetof_(T, field)))
// Restricted (non-aliased) pointers:
#define noalias __restrict__

// Placement new/delete are defined in <new> header for C++23
#include <new>

// Add raw data auto-casts to a structured or bitfield type
#define RawDataAccessors(structtype, rawtype) \
  structtype() { } \
  structtype(rawtype rawbits) { *((rawtype*)this) = rawbits; } \
  operator rawtype() const { return *((rawtype*)this); }

// Typecasts in bizarre ways required for binary form access
union W32orFloat { W32 w; float f; };
union W64orDouble {
  W64 w;
  double d;
  struct { W32 lo; W32s hi; } hilo;
  struct { W64 mantissa:52, exponent:11, negative:1; } ieee;
  // This format makes it easier to see if a NaN is a signalling NaN.
  struct { W64 mantissa:51, qnan:1, exponent:11, negative:1; } ieeenan;
};

static inline const float W32toFloat(W32 x) { union W32orFloat c; c.w = x; return c.f; }
static inline const W32 FloatToW32(float x) { union W32orFloat c; c.f = x; return c.w; }
static inline const double W64toDouble(W64 x) { union W64orDouble c; c.w = x; return c.d; }
static inline const W64 DoubleToW64(double x) { union W64orDouble c; c.d = x; return c.w; }

//
// Functional constructor
//

template <typename T> static inline T min(const T& a, const T& b) { return (a > b) ? b : a; }
template <typename T> static inline T max(const T& a, const T& b) { return (a > b) ? a : b; }
template <typename T> static inline T clipto(const T& v, const T& minv, const T& maxv) { return min(max(v, minv), maxv); }
// Use auto for return type deduction to handle mixed types
template <typename T, typename U, typename V>
static inline bool inrange(const T& v, const U& minv, const V& maxv) { return ((v >= minv) & (v <= maxv)); }
template <typename T> static inline T abs(T x) { return (x < 0) ? -x : x; }

// Bit fitting
static inline bool fits_in_signed_nbit(W64s v, int b) {
  return inrange(v, W64s(-(1ULL<< (b-1))), W64s(+(1ULL << (b-1))-1));
}

static inline bool fits_in_signed_nbit_tagged(W64s v, int b) {
  return inrange(v, W64s(-(1ULL<< (b-1))+1), W64s(+(1ULL << (b-1))-1));
}

static inline bool fits_in_signed_8bit(W64s v) { return fits_in_signed_nbit(v, 8); }
static inline bool fits_in_signed_16bit(W64s v) { return fits_in_signed_nbit(v, 16); }
static inline bool fits_in_signed_32bit(W64s v) { return fits_in_signed_nbit(v, 32); }

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define bit(x, n) (((x) >> (n)) & 1)

#define bitmask(l) (((l) == 64) ? (W64)(-1LL) : ((1LL << (l))-1LL))
#define bits(x, i, l) (((x) >> (i)) & bitmask(l))
#define lowbits(x, l) bits(x, 0, l)
#define setbit(x,i) ((x) |= (1LL << (i)))
#define clearbit(x, i) ((x) &= (W64)(~(1LL << (i))))
#define assignbit(x, i, v) ((x) = (((x) &= (W64)(~(1LL << (i)))) | (((W64)((bool)(v))) << i)));

#define foreach(i, n) for (size_t i = 0; i < (n); i++)

static inline W64s signext64(W64s x, const int i) { return (x << (64-i)) >> (64-i); }
static inline W32s signext32(W32s x, const int i) { return (x << (32-i)) >> (32-i); }
static inline W16s signext16(W16s x, const int i) { return (x << (16-i)) >> (16-i); }

static inline W64s bitsext64(W64s x, const int i, const int l) { return signext64(bits(x, i, l), l); }
static inline W32s bitsext32(W32s x, const int i, const int l) { return signext32(bits(x, i, l), l); }
static inline W16s bitsext16(W16s x, const int i, const int l) { return signext16(bits(x, i, l), l); }

// Vector types - platform independent definitions
// These are kept for API compatibility but simplified
typedef byte v16qi __attribute__ ((vector_size(16)));
typedef v16qi vec16b;
typedef W16 v8hi __attribute__ ((vector_size(16)));
typedef v8hi vec8w;
typedef float v4sf __attribute__ ((vector_size(16)));
typedef v4sf vec4f;
typedef W32 v4si __attribute__ ((vector_size(16)));
typedef v4si vec4i;
typedef float v2df __attribute__ ((vector_size(16)));
typedef v2df vec2d;

// Portable vector operations - use standard loops instead of SSE intrinsics
inline vec16b portable_pcmpeqb(vec16b a, vec16b b) {
  vec16b result;
  for (int i = 0; i < 16; i++) {
    result[i] = (a[i] == b[i]) ? 0xFF : 0x00;
  }
  return result;
}

inline vec8w portable_pcmpeqw(vec8w a, vec8w b) {
  vec8w result;
  for (int i = 0; i < 8; i++) {
    result[i] = (a[i] == b[i]) ? 0xFFFF : 0x0000;
  }
  return result;
}

inline vec4i portable_pcmpeqd(vec4i a, vec4i b) {
  vec4i result;
  for (int i = 0; i < 4; i++) {
    result[i] = (a[i] == b[i]) ? 0xFFFFFFFF : 0x00000000;
  }
  return result;
}

inline vec16b portable_psubusb(vec16b a, vec16b b) {
  vec16b result;
  for (int i = 0; i < 16; i++) {
    result[i] = (a[i] > b[i]) ? (a[i] - b[i]) : 0;
  }
  return result;
}

inline vec16b portable_paddusb(vec16b a, vec16b b) {
  vec16b result;
  for (int i = 0; i < 16; i++) {
    int sum = a[i] + b[i];
    result[i] = (sum > 255) ? 255 : sum;
  }
  return result;
}

inline vec16b portable_pandb(vec16b a, vec16b b) {
  vec16b result;
  for (int i = 0; i < 16; i++) {
    result[i] = a[i] & b[i];
  }
  return result;
}

inline W32 portable_pmovmskb(vec16b vec) {
  W32 mask = 0;
  for (int i = 0; i < 16; i++) {
    if (vec[i] & 0x80) mask |= (1 << i);
  }
  return mask;
}

// Compatibility aliases - point to portable implementations
#define x86_sse_pcmpeqb portable_pcmpeqb
#define x86_sse_pcmpeqw portable_pcmpeqw
#define x86_sse_pcmpeqd portable_pcmpeqd
#define x86_sse_psubusb portable_psubusb
#define x86_sse_paddusb portable_paddusb
#define x86_sse_pandb portable_pandb
#define x86_sse_pmovmskb portable_pmovmskb

inline vec16b x86_sse_psubusw(vec8w a, vec8w b) {
  vec8w result;
  for (int i = 0; i < 8; i++) {
    result[i] = (a[i] > b[i]) ? (a[i] - b[i]) : 0;
  }
  return (vec16b)result;
}

inline vec8w x86_sse_paddusw(vec8w a, vec8w b) {
  vec8w result;
  for (int i = 0; i < 8; i++) {
    int sum = a[i] + b[i];
    result[i] = (sum > 65535) ? 65535 : sum;
  }
  return result;
}

inline vec8w x86_sse_pandw(vec8w a, vec8w b) {
  vec8w result;
  for (int i = 0; i < 8; i++) {
    result[i] = a[i] & b[i];
  }
  return result;
}

inline vec16b x86_sse_packsswb(vec8w a, vec8w b) {
  vec16b result;
  for (int i = 0; i < 8; i++) {
    W16 val = a[i];
    if ((W16s)val < -128) result[i] = -128;
    else if ((W16s)val > 127) result[i] = 127;
    else result[i] = (byte)val;
  }
  for (int i = 0; i < 8; i++) {
    W16 val = b[i];
    if ((W16s)val < -128) result[i+8] = -128;
    else if ((W16s)val > 127) result[i+8] = 127;
    else result[i+8] = (byte)val;
  }
  return result;
}

inline W32 x86_sse_pmovmskw(vec8w vec) {
  return x86_sse_pmovmskb(x86_sse_packsswb(vec, vec)) & 0xff;
}

inline vec16b x86_sse_psadbw(vec16b a, vec16b b) {
  vec16b result = {};
  W64 sum0 = 0, sum1 = 0;
  for (int i = 0; i < 8; i++) {
    sum0 += std::abs(a[i] - b[i]);
  }
  for (int i = 8; i < 16; i++) {
    sum1 += std::abs(a[i] - b[i]);
  }
  // Store as 64-bit values in result
  W64* r = (W64*)&result;
  r[0] = sum0;
  r[1] = sum1;
  return result;
}

template <int i> inline W16 x86_sse_pextrw(vec16b a) {
  return ((W16*)&a)[i];
}

inline vec16b x86_sse_ldvbu(const vec16b* m) { return *m; }
inline void x86_sse_stvbu(vec16b* m, const vec16b ra) { *m = ra; }
inline vec8w x86_sse_ldvwu(const vec8w* m) { return *m; }
inline void x86_sse_stvwu(vec8w* m, const vec8w ra) { *m = ra; }

inline vec16b x86_sse_zerob() { vec16b rd = {}; return rd; }
inline vec16b x86_sse_onesb() { vec16b rd; memset(&rd, 0xFF, sizeof(rd)); return rd; }
inline vec8w x86_sse_zerow() { vec8w rd = {}; return rd; }
inline vec8w x86_sse_onesw() { vec8w rd; memset(&rd, 0xFF, sizeof(rd)); return rd; }

extern const byte byte_to_vec16b[256][16];
extern const byte index_bytes_vec16b[16][16];
extern const byte index_bytes_plus1_vec16b[16][16];

inline vec16b x86_sse_dupb(const byte b) {
  return *((vec16b*)&byte_to_vec16b[b]);
}

inline vec8w x86_sse_dupw(const W16 b) {
  W32 w = (b << 16) | b;
  vec8w v;
  W32* wp = (W32*)&v;
  wp[0] = w; wp[1] = w; wp[2] = w; wp[3] = w;
  return v;
}

// MXCSR handling - use portable fenv.h based approach
#include <cfenv>

inline void x86_set_mxcsr(W32 value) {
  // Map MXCSR rounding mode to fenv rounding mode
  int rc = (value >> 13) & 3;
  int femode;
  switch (rc) {
    case 0: femode = FE_TONEAREST; break;
    case 1: femode = FE_DOWNWARD; break;
    case 2: femode = FE_UPWARD; break;
    case 3: femode = FE_TOWARDZERO; break;
    default: femode = FE_TONEAREST; break;
  }
  std::fesetround(femode);
}

inline W32 x86_get_mxcsr() {
  int femode = std::fegetround();
  int rc;
  switch (femode) {
    case FE_TONEAREST: rc = 0; break;
    case FE_DOWNWARD: rc = 1; break;
    case FE_UPWARD: rc = 2; break;
    case FE_TOWARDZERO: rc = 3; break;
    default: rc = 0; break;
  }
  // Return default MXCSR with only rounding mode set
  return 0x1f80 | (rc << 13);
}

union MXCSR {
  struct { W32 ie:1, de:1, ze:1, oe:1, ue:1, pe:1, daz:1, im:1, dm:1, zm:1, om:1, um:1, pm:1, rc:2, fz:1; } fields;
  W32 data;

  MXCSR() { }
  MXCSR(W32 v) { data = v; }
  operator W32() const { return data; }
};
enum { MXCSR_ROUND_NEAREST, MXCSR_ROUND_DOWN, MXCSR_ROUND_UP, MXCSR_ROUND_TOWARDS_ZERO };
#define MXCSR_EXCEPTION_DISABLE_MASK 0x1f80 // OR this into mxcsr to disable all exceptions
#define MXCSR_DEFAULT 0x1f80 // default settings (no exceptions, defaults for rounding and denormals)

// Portable bit scan functions using C++20 <bit> header
inline W32 x86_bsf32(W32 b) { return (b == 0) ? 0 : std::countr_zero(b); }
inline W64 x86_bsf64(W64 b) { return (b == 0) ? 0 : std::countr_zero(b); }
inline W32 x86_bsr32(W32 b) { return (b == 0) ? 0 : 31 - std::countl_zero(b); }
inline W64 x86_bsr64(W64 b) { return (b == 0) ? 0 : 63 - std::countl_zero(b); }

// Portable bit test/set/reset/complement
template <typename T> inline bool x86_bt(T r, T b) { return (r >> b) & 1; }
template <typename T> inline bool x86_btn(T r, T b) { return !((r >> b) & 1); }

// Return the updated data; ignore the old value
template <typename T> inline W64 x86_bts(T r, T b) { return r | (T(1) << b); }
template <typename T> inline W64 x86_btr(T r, T b) { return r & ~(T(1) << b); }
template <typename T> inline W64 x86_btc(T r, T b) { return r ^ (T(1) << b); }

// Return the old value of the bit, but still update the data
template <typename T> inline bool x86_test_bts(T& r, T b) { bool c = (r >> b) & 1; r |= (T(1) << b); return c; }
template <typename T> inline bool x86_test_btr(T& r, T b) { bool c = (r >> b) & 1; r &= ~(T(1) << b); return c; }
template <typename T> inline bool x86_test_btc(T& r, T b) { bool c = (r >> b) & 1; r ^= (T(1) << b); return c; }

// Atomic bit operations using std::atomic
template <typename T> inline bool x86_locked_bts(T& r, T b) {
  std::atomic_ref<T> ar(r);
  T old = ar.fetch_or(T(1) << b);
  return (old >> b) & 1;
}
template <typename T> inline bool x86_locked_btr(T& r, T b) {
  std::atomic_ref<T> ar(r);
  T old = ar.fetch_and(~(T(1) << b));
  return (old >> b) & 1;
}
template <typename T> inline bool x86_locked_btc(T& r, T b) {
  std::atomic_ref<T> ar(r);
  T old = ar.fetch_xor(T(1) << b);
  return (old >> b) & 1;
}

// Portable byteswap using C++23
template <typename T> inline T bswap(T r) {
  if constexpr (sizeof(T) == 2) return std::byteswap(static_cast<uint16_t>(r));
  else if constexpr (sizeof(T) == 4) return std::byteswap(static_cast<uint32_t>(r));
  else if constexpr (sizeof(T) == 8) return std::byteswap(static_cast<uint64_t>(r));
  else return r;
}

static inline W16 x86_sse_maskeqb(const vec16b v, byte target) { return x86_sse_pmovmskb(x86_sse_pcmpeqb(v, x86_sse_dupb(target))); }

// This is a barrier for the compiler only, NOT the processor!
#define barrier() std::atomic_signal_fence(std::memory_order_seq_cst)

// Denote parallel sections for the compiler
#define parallel

// Portable atomic exchange
template <typename T>
static inline T xchg(T& v, T newv) {
  std::atomic_ref<T> av(v);
  return av.exchange(newv, std::memory_order_seq_cst);
}

// Portable atomic add
template <typename T>
static inline T xadd(T& v, T incr) {
  std::atomic_ref<T> av(v);
  return av.fetch_add(incr, std::memory_order_seq_cst);
}

// Portable compare-exchange
template <typename T>
static inline T cmpxchg(T& mem, T newv, T cmpv) {
  std::atomic_ref<T> amem(mem);
  amem.compare_exchange_strong(cmpv, newv, std::memory_order_seq_cst);
  // Return the old value in the slot (so we can check if it matches newv)
  return cmpv;
}

static inline void cpu_pause() {
#if defined(__x86_64__) || defined(__i386__)
  __builtin_ia32_pause();
#elif defined(__aarch64__)
  asm volatile("yield");
#else
  // Generic fallback - just a compiler barrier
  barrier();
#endif
}

static inline void prefetch(const void* x) {
  __builtin_prefetch(x, 0, 3);
}

// Portable CPUID - returns zeros on non-x86
static inline void cpuid(int op, W32& eax, W32& ebx, W32& ecx, W32& edx) {
#if defined(__x86_64__) || defined(__i386__)
  #if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__(
      "cpuid"
      : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
      : "a"(op), "c"(0)
    );
  #else
    eax = ebx = ecx = edx = 0;
  #endif
#else
  eax = ebx = ecx = edx = 0;
#endif
}

// Portable timestamp counter using chrono
static inline W64 rdtsc() {
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
}

// Portable rotate operations using C++20
template <typename T>
static inline T x86_ror(T r, int n) {
  return std::rotr(r, n);
}

template <typename T>
static inline T x86_rol(T r, int n) {
  return std::rotl(r, n);
}

template <typename T>
static inline T dupb(const byte b) { return T(b) * T(0x0101010101010101ULL); }

template <int n> struct lg { static const int value = 1 + lg<n/2>::value; };
template <> struct lg<1> { static const int value = 0; };
#define log2(v) (lg<(v)>::value)

template <int n> struct lg10 { static const int value = 1 + lg10<n/10>::value; };
template <> struct lg10<1> { static const int value = 0; };
template <> struct lg10<0> { static const int value = 0; };
#define log10(v) (lg10<(v)>::value)

template <int N, typename T>
static inline T foldbits(T a) {
  if (N == 0) return 0;

  const int B = (sizeof(T) * 8);
  const int S = (B / N) + ((B % N) ? 1 : 0);

  T z = 0;
  foreach (i, S) {
    z ^= a;
    a >>= N;
  }

  return lowbits(z, N);
}


// For specifying easy to read arrays
#define _ (0)

asmlinkage {
#include <unistd.h>
#include <sys/types.h>
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <fcntl.h>

#include <sys/mman.h>
#include <sys/utsname.h>
#include <sys/ptrace.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/user.h>
};

#include <cstdarg>

#include <syscalls.h>

#ifdef PAGE_SIZE
#undef PAGE_SIZE
#endif
// Standard page size - may vary by platform but 4096 is common
#define PAGE_SIZE 4096

/*
 * Make these math functions available even inside of member functions with the same name:
 */
static inline float fsqrt(float v) { return std::sqrt(v); }

template <typename T> static inline void setzero(T& x) { memset(&x, 0, sizeof(T)); }

#define HI32(x) (W32)((x) >> 32LL)
#define LO32(x) (W32)((x) & 0xffffffffLL)
#define CONCAT64(hi, lo) ((((W64)(hi)) << 32) + (((W64)(lo)) & 0xffffffffLL))

template <typename T, typename A> static inline T floor(T x, A a) { return (T)(((T)x) & ~((T)(a-1))); }
template <typename T, typename A> static inline T trunc(T x, A a) { return (T)(((T)x) & ~((T)(a-1))); }
template <typename T, typename A> static inline T ceil(T x, A a) { return (T)((((T)x) + ((T)(a-1))) & ~((T)(a-1))); }
template <typename T, typename A> static inline T mask(T x, A a) { return (T)(((T)x) & ((T)(a-1))); }

template <typename T, typename A> static inline T* floorptr(T* x, A a) { return (T*)floor((Waddr)x, a); }
template <typename T, typename A> static inline T* ceilptr(T* x, A a) { return (T*)ceil((Waddr)x, a); }
template <typename T, typename A> static inline T* maskptr(T* x, A a) { return (T*)mask((Waddr)x, a); }
static inline W64 mux64(W64 sel, W64 v0, W64 v1) { return (sel & v1) | ((~sel) & v0); }
template <typename T> static inline T mux(T sel, T v1, T v0) { return (sel & v1) | ((~sel) & v0); }

template <typename T> void swap(T& a, T& b) { T t = a;  a = b; b = t; }

//
// Portable branchless select using ternary operator
// The compiler will optimize this appropriately
//
template <typename T, typename K>
T select(K cond, T if0, T if1) {
  return cond ? if1 : if0;
}

template <typename T, typename K>
void condmove(K cond, T& v, T newv) {
  if (cond) v = newv;
}

#define typeof __typeof__
#define ptralign(ptr, bytes) ((typeof(ptr))((uintptr_t)(ptr) & ~((bytes)-1)))
#define ptrmask(ptr, bytes) ((typeof(ptr))((uintptr_t)(ptr) & ((bytes)-1)))

template <typename T>
inline void arraycopy(T* dest, const T* source, int count) { memcpy(dest, source, count * sizeof(T)); }

template <typename T, typename V>
inline void rawcopy(T& dest, const V& source) { memcpy(&dest, &source, sizeof(T)); }

static inline bool aligned(W64 address, int size) {
  return ((address & (W64)(size-1)) == 0);
}

inline bool strequal(const char* a, const char* b) {
  return (strcmp(a, b) == 0);
}

template <typename T, size_t size> size_t lengthof(T (&)[size]) { return size; }

extern const byte popcountlut8bit[];
extern const byte lsbindexlut8bit[];

static inline int popcount8bit(byte x) {
  return popcountlut8bit[x];
}

static inline int lsbindex8bit(byte x) {
  return lsbindexlut8bit[x];
}

// Use C++20 std::popcount
static inline int popcount(W32 x) {
  return std::popcount(x);
}

static inline int popcount64(W64 x) {
  return std::popcount(x);
}


extern const W64 expand_8bit_to_64bit_lut[256];

// LSB index using portable C++20
inline unsigned int lsbindex32(W32 n) { return (n == 0) ? 0 : std::countr_zero(n); }

inline int lsbindexi32(W32 n) {
  int r = lsbindex32(n);
  return (n ? r : -1);
}

inline unsigned int lsbindex64(W64 n) { return (n == 0) ? 0 : std::countr_zero(n); }

inline unsigned int lsbindexi64(W64 n) {
  int r = lsbindex64(n);
  return (n ? r : -1);
}

inline unsigned int lsbindex(W64 n) { return lsbindex64(n); }

// MSB index using portable C++20
inline unsigned int msbindex32(W32 n) { return (n == 0) ? 0 : 31 - std::countl_zero(n); }

inline int msbindexi32(W32 n) {
  int r = msbindex32(n);
  return (n ? r : -1);
}

inline unsigned int msbindex64(W64 n) { return (n == 0) ? 0 : 63 - std::countl_zero(n); }

inline unsigned int msbindexi64(W64 n) {
  int r = msbindex64(n);
  return (n ? r : -1);
}

inline unsigned int msbindex(W64 n) { return msbindex64(n); }

#define percent(x, total) (100.0 * ((float)(x)) / ((float)(total)))

inline int add_index_modulo(int index, int increment, int bufsize) {
  // Only if power of 2: return (index + increment) & (bufsize-1);
  index += increment;
  if (index < 0) index += bufsize;
  if (index >= bufsize) index -= bufsize;
  return index;
}

#include <superstl.h>

using namespace superstl;

ostream& operator <<(ostream& os, const vec16b& v);
ostream& operator ,(ostream& os, const vec16b& v);
ostream& operator <<(ostream& os, const vec8w& v);
ostream& operator ,(ostream& os, const vec8w& v);

#endif // __cplusplus

#endif // _GLOBALS_H_
