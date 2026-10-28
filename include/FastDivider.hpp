///
/// @file  FastDivider.hpp
/// @brief Fast division of 64-bit numerators < 2^63 by divisors
///        that are used many times. FastDivider replaces the
///        integer division by a multiplication with the divisor's
///        precomputed reciprocal and a right shift. This is much
///        faster than the CPU's integer division instruction,
///        especially on older CPUs.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FASTDIVIDER_HPP
#define FASTDIVIDER_HPP

#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>

#include <stdint.h>

#if defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_ARM64)) && \
    __has_include(<intrin.h>)
  #include <intrin.h>
#endif

namespace primecount {

/// Returns the upper 64 bits of the 128-bit product x * y
ALWAYS_INLINE uint64_t mulhi64(uint64_t x, uint64_t y)
{
#if defined(HAVE_INT128_T)
  return uint64_t((uint128_t(x) * y) >> 64);

#elif defined(_MSC_VER) && \
      (defined(_M_X64) || defined(_M_ARM64)) && \
      __has_include(<intrin.h>)
  return __umulh(x, y);

#else
  // Portable implementation for 32-bit CPUs, computes
  // x1 * y1 * 2^64 + (x1 * y0 + x0 * y1) * 2^32 + x0 * y0.
  uint32_t x0 = uint32_t(x);
  uint32_t x1 = uint32_t(x >> 32);
  uint32_t y0 = uint32_t(y);
  uint32_t y1 = uint32_t(y >> 32);
  uint64_t x0y0 = x0 * uint64_t(y0);
  uint64_t x0y1 = x0 * uint64_t(y1);
  uint64_t x1y0 = x1 * uint64_t(y0);
  uint64_t x1y1 = x1 * uint64_t(y1);
  uint64_t mid = x1y0 + (x0y0 >> 32);
  uint64_t mid_lo = uint32_t(mid) + x0y1;
  return x1y1 + (mid >> 32) + (mid_lo >> 32);
#endif
}

/// Computes (numhi * 2^64 + numlo) / den.
/// Requires numhi < den, so that the quotient fits into 64 bits.
///
inline uint64_t div128_64(uint64_t numhi,
                          uint64_t numlo,
                          uint64_t den)
{
  ASSERT(numhi < den);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))
  uint64_t quotient;
  uint64_t remainder;
  __asm__("div %[den]"
          : "=a"(quotient), "=d"(remainder)
          : [den] "r"(den), "a"(numlo), "d"(numhi));
  return quotient;

#elif defined(_MSC_VER) && \
      defined(_M_X64) && \
      _MSC_VER >= 1920 && \
      __has_include(<intrin.h>)
  uint64_t remainder;
  return _udiv128(numhi, numlo, den, &remainder);

#else
  // Knuth's Algorithm D in base 2^32, this is the same algorithm
  // as libdivide_128_div_64_to_64() in the libdivide library.
  // Normalize den so that its most significant bit is set.
  uint64_t shift = 63 - ilog2(den);
  den <<= shift;
  numhi <<= shift;
  if (shift > 0)
    numhi |= numlo >> (64 - shift);
  numlo <<= shift;

  uint64_t b = 1ull << 32;
  uint64_t den1 = den >> 32;
  uint64_t den0 = den & 0xffffffff;
  uint64_t num1 = numlo >> 32;
  uint64_t num0 = numlo & 0xffffffff;

  // Estimate the high quotient digit using den1,
  // the estimate is at most 2 too large.
  uint64_t qhat = numhi / den1;
  uint64_t rhat = numhi - qhat * den1;
  uint64_t c1 = qhat * den0;
  uint64_t c2 = rhat * b + num1;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint64_t q1 = uint32_t(qhat);

  // Estimate the low quotient digit using den1,
  // the estimate is at most 2 too large.
  uint64_t rem = numhi * b + num1 - q1 * den;
  qhat = rem / den1;
  rhat = rem - qhat * den1;
  c1 = qhat * den0;
  c2 = rhat * b + num0;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint64_t q0 = uint32_t(qhat);

  return (q1 << 32) | q0;
#endif
}

/// Divides 64-bit numerators n < 2^63 by a divisor d > 1 using
/// a multiplication by the divisor's precomputed reciprocal,
/// based on Ben Joffe's fast floor division algorithm:
/// https://www.benjoffe.com/fast-floordiv
///
/// For d > 1, s = floor(log2(d - 1)) and m = ceil(2^(64+s) / d)
/// the error e = m * d - 2^(64+s) < d <= 2^(s+1). Hence for
/// n < 2^63, n * e < 2^(64+s) and n / d = (n * m) >> (64 + s).
///
/// Packing FastDivider (9 bytes instead of 16) improves the cache
/// efficiency, this is faster even though some loads are unaligned.
///
#pragma pack(push, 1)
struct FastDivider
{
  FastDivider() = default;

  FastDivider(uint64_t d)
  {
    ASSERT(d > 1);
    shift = uint8_t(ilog2(d - 1));
    // multiplier = ceil(2^(64+shift) / d)
    multiplier = div128_64(1ull << shift, d - 1, d);
  }

  /// Largest supported numerator
  static constexpr uint64_t max()
  {
    return pstd::numeric_limits<int64_t>::max();
  }

  uint64_t multiplier;
  uint8_t shift;
};
#pragma pack(pop)

/// Requires n <= FastDivider::max()
ALWAYS_INLINE uint64_t operator/(uint64_t n, const FastDivider& d)
{
  ASSERT(n <= FastDivider::max());
  return mulhi64(n, d.multiplier) >> d.shift;
}

} // namespace

#endif
