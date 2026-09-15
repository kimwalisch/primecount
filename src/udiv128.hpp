///
/// @file  udiv128.hpp
/// @brief Integer division of small types is much faster than integer
///        division of large types on most CPUs. The fast_div(x, y)
///        function tries to take advantage of this by casting x and y
///        to smaller types (if possible) before doing the division.
///
///        On the x64 CPU architecture, if ENABLE_DIV32 is defined we
///        check at runtime if we can divide using the divl
///        instruction (64-bit / 32-bit = 32-bit) which is usually
///        faster than a full 64-bit division. On most CPUs before
///        2020 this significantly improves performance.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef UDIV128_HPP
#define UDIV128_HPP

#include <fast_div.hpp>
#include <ctz.hpp>
#include <macros.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <type_traits>

namespace primecount {

/// Used for (128-bit / 64-bit) = 64-bit.
/// This is the branchless correction variant of Knuth Algorithm D
/// used by libdivide's divllu() implementation, specialized to
/// a quotient that fits into 64 bits.
/// https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c
///
uint64_t udiv_128_by_64_to_64(uint64_t hi,
                              uint64_t lo,
                              uint64_t d)
{
  ASSERT(d != 0);
  ASSERT(hi < d);

  // 64-bit / 64-bit = 64-bit
  if (hi == 0)
    return lo / d;

  // 128-bit / 32-bit = 64-bit
  if (d <= UINT32_MAX)
  {
    uint64_t n1 = (hi << 32) | (lo >> 32);
    uint64_t n0 = lo;
    uint64_t dhi = d << 32;
    uint64_t q1 = n1 / d;
    uint64_t dividend = n0 - q1 * dhi;
    uint64_t q0 = dividend / d;

    return (q1 << 32) | q0;
  }

  // Normalize the divisor so that its most significant bit is set.
  // Since d >= 2^32, shift is in the range [0, 31].
  uint64_t shift = clz64(d);
  uint64_t den = d << shift;

  if (shift != 0)
  {
    hi = (hi << shift) | (lo >> (64 - shift));
    lo <<= shift;
  }

  // Work in base 2^32
  uint64_t den1 = den >> 32;
  uint64_t den0 = den & 0xffffffff;
  uint64_t n1 = lo >> 32;
  uint64_t n0 = lo & 0xffffffff;

  // Estimate and correct the high 32 quotient bits
  uint64_t q1 = hi / den1;
  uint64_t rhat = hi % den1;
  uint64_t c1 = q1 * den0;
  uint64_t c2 = (rhat << 32) + n1;

  if (c1 > c2)
    q1 -= (c1 - c2 > den) ? 2 : 1;

  // Compute the true partial remainder
  uint64_t rem = (hi << 32) + n1 - q1 * den;

  // Estimate and correct the low 32 quotient bits
  uint64_t q0 = rem / den1;
  rhat = rem % den1;
  c1 = q0 * den0;
  c2 = (rhat << 32) + n0;

  if (c1 > c2)
    q0 -= (c1 - c2 > den) ? 2 : 1;

  return (q1 << 32) | q0;
}

/// Used for (128-bit / 64-bit) = 128-bit.
/// Handles general 128-bit division and uses optimized narrowing
/// division when the divisor and quotient fit into 64 bits. The general
/// 128/64-bit narrowing case uses the branchless correction variant
/// of Knuth Algorithm D from libdivide's divllu() implementation.
/// https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c
///
uint128_t udiv_128_by_64_to_128(uint64_t hi,
                                uint64_t lo,
                                uint64_t d)
{
  // Prevent division by zero
  ASSERT(d != 0);

  // 128-bit / 64-bit = 128-bit
  if (hi >= d)
  {
    uint128_t n = (uint128_t(hi) << 64) | lo;
    return n / d;
  }

  // 64-bit / 64-bit = 64-bit
  if (hi == 0)
    return lo / d;

  // 128-bit / 32-bit = 64-bit
  if (d <= UINT32_MAX)
  {
    uint64_t n1 = (hi << 32) | (lo >> 32);
    uint64_t n0 = lo;
    uint64_t dhi = d << 32;
    uint64_t q1 = n1 / d;
    uint64_t dividend = n0 - q1 * dhi;
    uint64_t q0 = dividend / d;

    return (q1 << 32) | q0;
  }

  // Normalize the divisor so that its most significant bit is set.
  // Since d >= 2^32, shift is in the range [0, 31].
  uint64_t shift = clz64(d);
  uint64_t den = d << shift;

  if (shift != 0)
  {
    hi = (hi << shift) | (lo >> (64 - shift));
    lo <<= shift;
  }

  // Work in base 2^32
  uint64_t den1 = den >> 32;
  uint64_t den0 = den & 0xffffffff;
  uint64_t n1 = lo >> 32;
  uint64_t n0 = lo & 0xffffffff;

  // Estimate and correct the high 32 quotient bits
  uint64_t q1 = hi / den1;
  uint64_t rhat = hi % den1;
  uint64_t c1 = q1 * den0;
  uint64_t c2 = (rhat << 32) + n1;

  if (c1 > c2)
    q1 -= (c1 - c2 > den) ? 2 : 1;

  // Compute the true partial remainder
  uint64_t rem = (hi << 32) + n1 - q1 * den;

  // Estimate and correct the low 32 quotient bits
  uint64_t q0 = rem / den1;
  rhat = rem % den1;
  c1 = q0 * den0;
  c2 = (rhat << 32) + n0;

  if (c1 > c2)
    q0 -= (c1 - c2 > den) ? 2 : 1;

  return (q1 << 32) | q0;
}

/// Used for (128-bit / 128-bit) = 128-bit.
/// Handles general 128-bit division and uses optimized narrowing
/// division when the divisor and quotient fit into 64 bits. The general
/// 128/64-bit narrowing case uses the branchless correction variant
/// of Knuth Algorithm D from libdivide's divllu() implementation.
/// https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c
///
uint128_t udiv_128_by_128_to_128(uint64_t hi,
                                 uint64_t lo,
                                 uint64_t d1,
                                 uint64_t d0)
{
  // Prevent division by zero
  ASSERT(d0 != 0 || d1 != 0);

  // 128-bit / 128-bit = 128-bit
  if (d1 != 0 || hi >= d0)
  {
    uint128_t n = (uint128_t(hi) << 64) | lo;
    uint128_t d = (uint128_t(d1) << 64) | d0;
    return n / d;
  }

  // 64-bit / 64-bit = 64-bit
  if (hi == 0)
    return lo / d0;

  // 128-bit / 32-bit = 64-bit
  if (d0 <= UINT32_MAX)
  {
    uint64_t n1 = (hi << 32) | (lo >> 32);
    uint64_t n0 = lo;
    uint64_t dhi = d0 << 32;
    uint64_t q1 = n1 / d0;
    uint64_t dividend = n0 - q1 * dhi;
    uint64_t q0 = dividend / d0;

    return (q1 << 32) | q0;
  }

  // Normalize the divisor so that its most significant bit is set.
  // Since d0 >= 2^32, shift is in the range [0, 31].
  uint64_t shift = clz64(d0);
  uint64_t den = d0 << shift;

  if (shift != 0)
  {
    hi = (hi << shift) | (lo >> (64 - shift));
    lo <<= shift;
  }

  // Work in base 2^32
  uint64_t den1 = den >> 32;
  uint64_t den0 = den & 0xffffffff;
  uint64_t n1 = lo >> 32;
  uint64_t n0 = lo & 0xffffffff;

  // Estimate and correct the high 32 quotient bits
  uint64_t q1 = hi / den1;
  uint64_t rhat = hi % den1;
  uint64_t c1 = q1 * den0;
  uint64_t c2 = (rhat << 32) + n1;

  if (c1 > c2)
    q1 -= (c1 - c2 > den) ? 2 : 1;

  // Compute the true partial remainder
  uint64_t rem = (hi << 32) + n1 - q1 * den;

  // Estimate and correct the low 32 quotient bits
  uint64_t q0 = rem / den1;
  rhat = rem % den1;
  c1 = q0 * den0;
  c2 = (rhat << 32) + n0;

  if (c1 > c2)
    q0 -= (c1 - c2 > den) ? 2 : 1;

  return (q1 << 32) | q0;
}

} // namespace

#endif
