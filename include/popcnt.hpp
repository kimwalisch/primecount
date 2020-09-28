///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits inside
///        an array or a 64-bit word.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <stdint.h>

#if defined(__has_include)
  #define HAS_INCLUDE(header) __has_include(header)
#else
  // If the __has_include() macro does not exist
  // we assume that the header file exists.
  #define HAS_INCLUDE(header) 1
#endif

#if defined(ENABLE_POPCNT)

#if !defined(__has_builtin)
  #define __has_builtin(x) 0
#endif

// GCC & Clang
#if defined(__GNUC__) || \
    __has_builtin(__builtin_popcountll)

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  return __builtin_popcountll(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_WIN64) && \
      HAS_INCLUDE(<nmmintrin.h>)

#include <nmmintrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  return _mm_popcnt_u64(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_WIN32) && \
      HAS_INCLUDE(<nmmintrin.h>)

#include <nmmintrin.h>

namespace {

inline uint64_t popcnt64(uint64_t x)
{
  return _mm_popcnt_u32((uint32_t) x) +
         _mm_popcnt_u32((uint32_t)(x >> 32));
}

} // namespace

#else

// Since hardware popcount support is very important for
// performance in primecount we want the compilation to fail in
// case the compiler does not select any of the implementations
// above (which use compiler intrinsics to enable hardware
// popcount support). This way we can be relatively sure that
// when primecount has been compiled successfully, primecount
// will have hardware popcount support.
#error "No fast popcount function implemented for this compiler!"

#endif
#endif

namespace {

#if !defined(ENABLE_POPCNT)

/// This uses fewer arithmetic operations than any other known
/// implementation on machines with fast multiplication.
/// It uses 12 arithmetic operations, one of which is a multiply.
/// https://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
///
inline uint64_t popcnt64(uint64_t x)
{
  uint64_t m1 = 0x5555555555555555ull;
  uint64_t m2 = 0x3333333333333333ull;
  uint64_t m4 = 0x0F0F0F0F0F0F0F0Full;
  uint64_t h01 = 0x0101010101010101ull;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

#endif

} // namespace

#endif // POPCNT_HPP
