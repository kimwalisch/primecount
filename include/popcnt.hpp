///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits inside
///        an array or a 64-bit word.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <CPUID.hpp>
#include <stdint.h>

#if !defined(__has_builtin)
  #define __has_builtin(x) 0
#endif

#if !defined(__has_include)
  #define __has_include(x) 0
#endif

namespace {

/// This uses fewer arithmetic operations than any other known
/// implementation on machines with fast multiplication.
/// It uses 12 arithmetic operations, one of which is a multiply.
/// http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
///
inline uint64_t popcnt64_bitwise(uint64_t x)
{
  uint64_t m1 = 0x5555555555555555ll;
  uint64_t m2 = 0x3333333333333333ll;
  uint64_t m4 = 0x0F0F0F0F0F0F0F0Fll;
  uint64_t h01 = 0x0101010101010101ll;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

} // namespace

// GCC & Clang
#if defined(__GNUC__) || \
    __has_builtin(__builtin_popcountl)

// CPUID is only enabled on x86 and x86-64 CPUs
// if the user compiles without -mpopcnt.
#if defined(ENABLE_CPUID_POPCNT)

#if defined(__x86_64__)

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  // On my AMD EPYC 7642 CPU using GCC 12 this runtime
  // check incurs an overall overhead of about 1%.
  if_likely(HAS_CPUID_POPCNT)
  {
    __asm__("popcnt %1, %0" : "=r"(x) : "r"(x));
    return x;
  }
  else
  {
    // On x86 and x64 CPUs when using the GCC compiler
    // __builtin_popcount*(x) is slow (not inlined function call)
    // when compiling without -mpopcnt. Therefore we avoid
    // using __builtin_popcount*(x) here.
    return popcnt64_bitwise(x);
  }
}

} // namespace

#elif defined(__i386__)

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  if_likely(HAS_CPUID_POPCNT)
  {
    uint32_t x0 = uint32_t(x);
    uint32_t x1 = uint32_t(x >> 32);
    __asm__("popcnt %1, %0" : "=r"(x0) : "r"(x0));
    __asm__("popcnt %1, %0" : "=r"(x1) : "r"(x1));
    return x0 + x1;
  }
  else
  {
    // On x86 and x64 CPUs when using the GCC compiler
    // __builtin_popcount*(x) is slow (not inlined function call)
    // when compiling without -mpopcnt. Therefore we avoid
    // using __builtin_popcount*(x) here.
    return popcnt64_bitwise(x);
  }
}

} // namespace

#endif // i386

#else // GCC & Clang (no CPUID, not x86)

namespace {

inline uint64_t popcnt64(uint64_t x)
{
#if __cplusplus >= 201703L
  if constexpr(sizeof(int) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcount(x);
  else if constexpr(sizeof(long) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcountl(x);
  else if constexpr(sizeof(long long) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcountll(x);
#else
    return (uint64_t) __builtin_popcountll(x);
#endif
}

} // namespace

#endif // GCC & Clang

#elif defined(_MSC_VER) && \
      defined(_M_X64) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
#if defined(HAS_POPCNT)
  return __popcnt64(x);
#elif defined(ENABLE_CPUID_POPCNT)
  if_likely(HAS_CPUID_POPCNT)
    return __popcnt64(x);
  else
    return popcnt64_bitwise(x);
#else
  return popcnt64_bitwise(x);
#endif
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_M_IX86) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
#if defined(HAS_POPCNT)
  return __popcnt(uint32_t(x)) +
         __popcnt(uint32_t(x >> 32));
#elif defined(ENABLE_CPUID_POPCNT)
  if_likely(HAS_CPUID_POPCNT)
    return __popcnt(uint32_t(x)) +
           __popcnt(uint32_t(x >> 32));
  else
    return popcnt64_bitwise(x);
#else
  return popcnt64_bitwise(x);
#endif
}

} // namespace

#elif __cplusplus >= 202002L

#include <bit>

namespace {

/// We only use the C++ standard library as a fallback if there
/// are no compiler intrinsics available for POPCNT.
/// Compiler intrinsics often generate faster assembly.
inline uint64_t popcnt64(uint64_t x)
{
  return std::popcount(x);
}

} // namespace

#endif

#endif // POPCNT_HPP
