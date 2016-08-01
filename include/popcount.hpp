///
/// @file    popcount.hpp
/// @brief   Functions to count the number of 1 bits inside a 64-bit
///          word or an array. If HAVE_POPCNT is defined then
///          popcnt_u64(x) will be used which uses the POPCNT
///          instruction (requires SSE4.2 for x86). For performance
///          reasons all algorithms are defined inline.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCOUNT_HPP
#define POPCOUNT_HPP

#if !defined(__STDC_CONSTANT_MACROS)
  #define __STDC_CONSTANT_MACROS
#endif

#include <stdint.h>

// HAVE_POPCNT is defined if the CPU supports the POPCNT instruction
#if defined(HAVE_POPCNT)

#if defined(_MSC_VER) && defined(_WIN64)
#define HAVE_POPCNT_U64

#include <nmmintrin.h>

namespace primecount {

inline uint64_t popcnt_u64(uint64_t x)
{
  return _mm_popcnt_u64(x);
}

}

#elif defined(_MSC_VER) && defined(_WIN32)
#define HAVE_POPCNT_U64

#include <nmmintrin.h>

namespace primecount {

inline uint64_t popcnt_u64(uint64_t x)
{
  return _mm_popcnt_u32((uint32_t) x) + 
         _mm_popcnt_u32((uint32_t)(x >> 32));
}

}

#elif defined(HAVE___BUILTIN_POPCOUNT) && defined(__i386__)
#define HAVE_POPCNT_U64

namespace primecount {

inline uint64_t popcnt_u64(uint64_t x)
{
  return __builtin_popcount((uint32_t) x) +
         __builtin_popcount((uint32_t)(x >> 32));
}

}

#elif defined(HAVE___BUILTIN_POPCOUNTLL)
#define HAVE_POPCNT_U64

namespace primecount {

inline uint64_t popcnt_u64(uint64_t x)
{
  // currently the type `long long' is 64-bits on all available CPUs 
  return __builtin_popcountll(x);
}

}

#else /* Error */

#if defined(_MSC_VER)
  #pragma error( "POPCNT not supported, remove \"/D HAVE_POPCNT\" from Makefile.msvc" )
#else
  #error "POPCNT not supported, use --disable-popcnt"
#endif

#endif
#endif /* HAVE_POPCNT */

#if defined(HAVE_POPCNT_U64)

namespace primecount {

inline uint64_t popcount_u64(uint64_t x)
{
  return popcnt_u64(x);
}

/// Count the number of 1 bits in an array using the POPCNT
/// instruction. On x86 CPUs this requires SSE4.2.
///
inline uint64_t popcount_u64(const uint64_t* array, uint64_t size)
{
  uint64_t sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
  uint64_t limit = size - size % 4;
  uint64_t i = 0;

  for (; i < limit; i += 4)
  {
    sum0 += popcnt_u64(array[i+0]);
    sum1 += popcnt_u64(array[i+1]);
    sum2 += popcnt_u64(array[i+2]);
    sum3 += popcnt_u64(array[i+3]);
  }

  uint64_t total = sum0 + sum1 + sum2 + sum3;

  for (; i < size; i++)
    total += popcnt_u64(array[i]);

  return total;
}

} // namespace

#else /* no POPCNT */

namespace primecount {

/// This uses fewer arithmetic operations than any other known  
/// implementation on machines with fast multiplication.
/// It uses 12 arithmetic operations, one of which is a multiply.
/// http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
///
inline uint64_t popcount_u64(uint64_t x)
{
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t h01 = UINT64_C(0x0101010101010101);

  x -=            (x >> 1)  & m1;
  x = (x & m2) + ((x >> 2)  & m2);
  x = (x +        (x >> 4)) & m4;
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
inline uint64_t popcount_u64(const uint64_t* array, uint64_t size)
{
  uint64_t total = 0;
  uint64_t ones = 0, twos = 0, fours = 0, eights = 0;
  uint64_t twosA, twosB, foursA, foursB;
  uint64_t limit = size - size % 8;
  uint64_t i = 0;

  for(; i < limit; i += 8)
  {
    CSA(twosA, ones, ones, array[i+0], array[i+1]);
    CSA(twosB, ones, ones, array[i+2], array[i+3]);
    CSA(foursA, twos, twos, twosA, twosB);    
    CSA(twosA, ones, ones, array[i+4], array[i+5]);
    CSA(twosB, ones, ones, array[i+6], array[i+7]);
    CSA(foursB, twos, twos, twosA, twosB);
    CSA(eights, fours, fours, foursA, foursB);

    total += popcount_u64(eights);
  }

  total *= 8;
  total += 4 * popcount_u64(fours);
  total += 2 * popcount_u64(twos);
  total += 1 * popcount_u64(ones);

  for(; i < size; i++)
    total += popcount_u64(array[i]);

  return total;
}

} // namespace

#endif

#endif /* POPCOUNT_HPP */
