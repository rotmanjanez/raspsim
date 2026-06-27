#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <cstdint>

namespace x86sim {

typedef __SIZE_TYPE__ size_t;
// W64 is std::uint64_t so it is the same type as word_t (the register-file
// element type). On LP64 that is unsigned long, on LLP64 unsigned long long;
// spelling it out as either breaks reference binding on the other ABI.
typedef std::uint64_t W64;
typedef std::int64_t W64s;
typedef unsigned int W32;
typedef signed int W32s;
typedef unsigned short W16;
typedef signed short W16s;
typedef unsigned char byte;
typedef unsigned char W8;
typedef signed char W8s;
#define null NULL

typedef W64 Waddr;


} // namespace x86sim

#endif