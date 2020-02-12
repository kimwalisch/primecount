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

#if !defined(ENABLE_POPCNT)

namespace {

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

} // namespace

#elif (defined(__ARM_NEON) || \
       defined(__aarch64__)) && \
       HAS_INCLUDE(<arm_neon.h>)

#include <algorithm>
#include <arm_neon.h>

namespace {

inline uint64x2_t vpadalq(uint64x2_t sum, uint8x16_t t)
{
  return vpadalq_u32(sum, vpaddlq_u16(vpaddlq_u8(t)));
}

/// Count the number of 1 bits in the data array.
/// ARM NEON has a vector popcount instruction but no scalar
/// popcount instruction that's why the ARM popcount function
/// is more complicated than the x86 popcount function.
/// This function has been copied from:
/// https://github.com/kimwalisch/libpopcnt
///
inline uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t i = 0;
  uint64_t cnt = 0;

  if (size >= 8)
  {
    const uint8_t* data8 = (const uint8_t*) data;
    uint64_t bytes = size * sizeof(uint64_t);
    uint64_t chunk_size = 64;
    uint64_t iters = bytes / chunk_size;

    uint64x2_t sum = vcombine_u64(vcreate_u64(0), vcreate_u64(0));
    uint8x16_t zero = vcombine_u8(vcreate_u8(0), vcreate_u8(0));

    do
    {
      uint8x16_t t0 = zero;
      uint8x16_t t1 = zero;
      uint8x16_t t2 = zero;
      uint8x16_t t3 = zero;

      // After every 31 iterations we need to add the
      // temporary sums (t0, t1, t2, t3) to the total sum.
      // We must ensure that the temporary sums <= 255
      // and 31 * 8 bits = 248 which is OK.
      uint64_t limit = std::min(i + 31, iters);

      // Each iteration processes 64 bytes
      for (; i < limit; i++)
      {
        uint8x16x4_t input = vld4q_u8(data8);
        data8 += chunk_size;

        t0 = vaddq_u8(t0, vcntq_u8(input.val[0]));
        t1 = vaddq_u8(t1, vcntq_u8(input.val[1]));
        t2 = vaddq_u8(t2, vcntq_u8(input.val[2]));
        t3 = vaddq_u8(t3, vcntq_u8(input.val[3]));
      }

      sum = vpadalq(sum, t0);
      sum = vpadalq(sum, t1);
      sum = vpadalq(sum, t2);
      sum = vpadalq(sum, t3);
    }
    while (i < iters);

    uint64_t tmp[2];
    vst1q_u64(tmp, sum);
    cnt += tmp[0];
    cnt += tmp[1];

    uint64_t bytes_processed = iters * chunk_size;
    i = bytes_processed / sizeof(uint64_t);
  }

  // Process the remaining bytes
  for (; i < size; i++)
    cnt += popcnt64(data[i]);

  return cnt;
}

} // namespace

#else

namespace {

/// Used for all CPUs except ARM
inline uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t cnt = 0;

  for (uint64_t i = 0; i < size; i++)
    cnt += popcnt64(data[i]);

  return cnt;
}

} // namespace

#endif

#endif // POPCNT_HPP
