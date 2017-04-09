///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits inside
///        an array or a 64-bit word.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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

#if defined(__GNUC__) || \
    __has_builtin(__builtin_popcountll)

inline uint64_t popcnt64(uint64_t x)
{
  return __builtin_popcountll(x);
}

#elif defined(_MSC_VER) && \
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

#else

inline uint64_t popcnt64(uint64_t x)
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

#endif

inline uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t cnt = 0;
  uint64_t i = 0;
  uint64_t limit = size - size % 4;

  for (; i < limit; i += 4)
  {
    cnt += popcnt64(data[i+0]);
    cnt += popcnt64(data[i+1]);
    cnt += popcnt64(data[i+2]);
    cnt += popcnt64(data[i+3]);
  }

  for (; i < size; i++)
    cnt += popcnt64(data[i]);

  return cnt;
}

#endif /* POPCNT_HPP */
