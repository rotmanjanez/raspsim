#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <cstdint>
#include <cstddef>

// Fixed-width integer types using standard C++ types
typedef std::size_t size_t;
typedef uint64_t W64;
typedef int64_t W64s;
typedef uint32_t W32;
typedef int32_t W32s;
typedef uint16_t W16;
typedef int16_t W16s;
typedef uint8_t byte;
typedef uint8_t W8;
typedef int8_t W8s;

#define null nullptr

// Platform-independent address type
#if INTPTR_MAX == INT64_MAX
typedef W64 Waddr;
#else
typedef W32 Waddr;
#endif

#endif
