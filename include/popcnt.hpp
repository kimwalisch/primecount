///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits
///        inside a 64-bit word.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <stdint.h>

#ifndef __has_builtin
  #define __has_builtin(x) 0
#endif

#if defined(_MSC_VER) && \
    defined(_WIN64)

#include <nmmintrin.h>

inline uint64_t popcnt64(uint64_t x)
{
  return _mm_popcnt_u64(x);
}

#elif defined(_MSC_VER) && \
      defined(_WIN32)

#include <nmmintrin.h>

inline uint64_t popcnt64(uint64_t x)
{
  return _mm_popcnt_u32((uint32_t) x) + 
         _mm_popcnt_u32((uint32_t)(x >> 32));
}

#elif defined(__GNUC__) || \
      __has_builtin(__builtin_popcountll)

inline uint64_t popcnt64(uint64_t x)
{
  return __builtin_popcountll(x);
}

#else

inline uint64_t popcnt64(uint64_t x)
{
  const uint64_t m1 = 0x5555555555555555ll;
  const uint64_t m2 = 0x3333333333333333ll;
  const uint64_t m4 = 0x0F0F0F0F0F0F0F0Fll;
  const uint64_t h01 = 0x0101010101010101ll;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

#endif

#endif /* POPCNT_HPP */
