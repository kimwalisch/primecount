///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits inside
///        an array or a 64-bit word.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <stdint.h>

#if !defined(DISABLE_POPCNT)

#ifndef __has_builtin
  #define __has_builtin(x) 0
#endif

#if defined(__GNUC__) || __has_builtin(__builtin_popcountll)

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

// fallback mode
#define DISABLE_POPCNT

#endif
#endif

#if !defined(DISABLE_POPCNT)

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

#endif

#if defined(DISABLE_POPCNT)

/// This uses fewer arithmetic operations than any other known
/// implementation on machines with fast multiplication.
/// It uses 12 arithmetic operations, one of which is a multiply.
/// https://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
///
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

/// Carry-save adder (CSA).
/// @see Chapter 5 in "Hacker's Delight" 2nd edition.
///
inline void CSA(uint64_t& h, uint64_t& l, uint64_t a, uint64_t b, uint64_t c)
{
  uint64_t u = a ^ b;
  h = (a & b) | (u & c);
  l = u ^ c;
}

/// Harley-Seal popcount (3rd iteration).
/// The Harley-Seal popcount algorithm is one of the fastest algorithms
/// for counting 1 bits in an array using only integer operations.
/// This implementation uses only 6.38 instructions per 64-bit word.
/// @see Chapter 5 in "Hacker's Delight" 2nd edition.
///
inline uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t cnt = 0;
  uint64_t ones = 0, twos = 0, fours = 0, eights = 0;
  uint64_t twosA, twosB, foursA, foursB;
  uint64_t limit = size - size % 8;
  uint64_t i = 0;

  for(; i < limit; i += 8)
  {
    CSA(twosA, ones, ones, data[i+0], data[i+1]);
    CSA(twosB, ones, ones, data[i+2], data[i+3]);
    CSA(foursA, twos, twos, twosA, twosB);
    CSA(twosA, ones, ones, data[i+4], data[i+5]);
    CSA(twosB, ones, ones, data[i+6], data[i+7]);
    CSA(foursB, twos, twos, twosA, twosB);
    CSA(eights, fours, fours, foursA, foursB);

    cnt += popcnt64(eights);
  }

  cnt *= 8;
  cnt += 4 * popcnt64(fours);
  cnt += 2 * popcnt64(twos);
  cnt += 1 * popcnt64(ones);

  for(; i < size; i++)
    cnt += popcnt64(data[i]);

  return cnt;
}

#endif

#endif // POPCNT_HPP
