///
/// @file  FastDivider.hpp
/// @brief Fast division of 64-bit numerators < 2^63 by divisors
///        that are used many times. FastDivider replaces the
///        integer division by a multiplication with the divisor's
///        precomputed reciprocal and a right shift. This is much
///        faster than the CPU's integer division instruction,
///        especially on older CPUs. FastDivider is not supported
///        on 32-bit CPUs, cmake/FastDivider.cmake checks whether
///        the compiler supports FastDivider.
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
#include <type_traits>

#if defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_ARM64)) && \
    __has_include(<intrin.h>)
  #include <intrin.h>
#endif

namespace primecount {

/// Returns the upper 64 bits of the 128-bit product x * y
ALWAYS_INLINE uint64_t mulhi_u64(uint64_t x, uint64_t y)
{
#if defined(HAVE_INT128_T)
  return uint64_t((uint128_t(x) * y) >> 64);

#elif defined(_MSC_VER) && \
      (defined(_M_X64) || defined(_M_ARM64)) && \
      __has_include(<intrin.h>)
  return __umulh(x, y);

#else
  #error "FastDivider.hpp requires int128_t or MSVC's __umulh()"
#endif
}

/// Computes (numhi * 2^64 + numlo) / den.
/// Requires numhi < den, so that the quotient fits into 64 bits.
///
inline uint64_t udiv128_64(uint64_t numhi,
                           uint64_t numlo,
                           uint64_t den)
{
  ASSERT(numhi < den);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || \
    defined(__clang__))
  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(numlo), "+d"(numhi) : [divider] "r"(den));
  return numlo;

#elif defined(HAVE_INT128_T)
  uint128_t numer = (uint128_t(numhi) << 64) | numlo;
  return uint64_t(numer / den);

#else
  // This code has been copied and adapted from libdivide's
  // divlu.c file which is released into the public domain:
  // https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c

  // We work in base 2**32.
  // A uint32 holds a single digit. A uint64 holds two digits.
  // Our numerator is conceptually [num3, num2, num1, num0].
  // Our denominator is [den1, den0].
  const uint64_t b = 1ull << 32;

  // Determine the normalization factor. We multiply den by this,
  // so that its leading digit is at least half b. In binary this
  // means just shifting left by the number of leading zeros, so
  // that there's a 1 in the MSB.
  // We also shift numer by the same amount. This cannot overflow
  // because numhi < den. The expression (-shift & 63) is the same
  // as (64 - shift), except it avoids the UB of shifting by 64.
  // The funny bitwise 'and' ensures that numlo does not get shifted
  // into numhi if shift is 0. clang 11 has an x86 codegen bug here:
  // see LLVM bug 50118. The sequence below avoids it.
  int shift = int(63 - ilog2(den));
  den <<= shift;
  numhi <<= shift;
  numhi |= (numlo >> (-shift & 63)) & uint64_t(-int64_t(shift) >> 63);
  numlo <<= shift;

  // Extract the low digits of the numerator and
  // both digits of the denominator.
  uint32_t num1 = uint32_t(numlo >> 32);
  uint32_t num0 = uint32_t(numlo & 0xffffffff);
  uint32_t den1 = uint32_t(den >> 32);
  uint32_t den0 = uint32_t(den & 0xffffffff);

  // We wish to compute q1 = [n3 n2 n1] / [d1 d0].
  // Estimate q1 as [n3 n2] / [d1], and then correct it.
  // Note while qhat may be 2 digits, q1 is always 1 digit.
  uint64_t qhat = numhi / den1;
  uint64_t rhat = numhi % den1;
  uint64_t c1 = qhat * den0;
  uint64_t c2 = rhat * b + num1;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint32_t q1 = uint32_t(qhat);

  // Compute the true (partial) remainder
  uint64_t rem = numhi * b + num1 - q1 * den;

  // We wish to compute q0 = [rem1 rem0 n0] / [d1 d0].
  // Estimate q0 as [rem1 rem0] / [d1] and correct it.
  qhat = rem / den1;
  rhat = rem % den1;
  c1 = qhat * den0;
  c2 = rhat * b + num0;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint32_t q0 = uint32_t(qhat);

  return (uint64_t(q1) << 32) | q0;
#endif
}

/// Divides 64-bit numerators n < 2^63 by a divisor d > 1 using
/// a multiplication by the divisor's precomputed reciprocal,
/// based on Granlund and Montgomery's "Division by Invariant
/// Integers using Multiplication" (1994):
/// https://gmplib.org/~tege/divcnst-pldi94.pdf
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
    multiplier = udiv128_64(1ull << shift, d - 1, d);
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
template <typename T>
ALWAYS_INLINE typename std::enable_if<pstd::is_same<T, uint64_t>::value, T>::type
operator/(T n, const FastDivider& d)
{
  ASSERT(n <= FastDivider::max());
  return mulhi_u64(n, d.multiplier) >> d.shift;
}

} // namespace

#endif
