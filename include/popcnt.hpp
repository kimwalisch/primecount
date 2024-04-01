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

#include <CpuID.hpp>
#include <stdint.h>

#if !defined(__has_builtin)
  #define __has_builtin(x) 0
#endif

#if !defined(__has_include)
  #define __has_include(x) 0
#endif

// GCC & Clang
#if defined(__GNUC__) || \
    __has_builtin(__builtin_popcountl)

namespace {

inline uint64_t popcnt64(uint64_t x)
{
// CPUID is only enabled on x86 and x86-64 CPUs
// if the user compiles without -mpopcnt.
#if defined(HAS_CPUID_POPCNT)
#if defined(__x86_64__)
  // On my AMD EPYC 7642 CPU using GCC 12 this runtime
  // check incurs an overall overhead of 2.5%.
  if_likely(CPUID_POPCNT)
  {
    __asm__("popcnt %1, %0" : "=r"(x) : "r"(x));
    return x;
  }
#elif defined(__i386__)
  if_likely(CPUID_POPCNT)
  {
    uint32_t x0 = uint32_t(x);
    uint32_t x1 = uint32_t(x >> 32);
    __asm__("popcnt %1, %0" : "=r"(x0) : "r"(x0));
    __asm__("popcnt %1, %0" : "=r"(x1) : "r"(x1));
    return x0 + x1;
  }
#endif
#endif

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

#elif defined(_MSC_VER) && \
      defined(_M_X64) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  return __popcnt64(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_M_IX86) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  return __popcnt((uint32_t) x) +
         __popcnt((uint32_t)(x >> 32));
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
