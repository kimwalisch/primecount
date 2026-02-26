///
/// @file  fast_div_128_64.hpp
/// @brief This file contains the function fast_div_128_64() that
///        is a modified version of the libdivide_128_div_64_to_64()
///        function from the libdivide library,
///        https://github.com/ridiculousfish/libdivide.
///        Hence this file is licensed under the same zlib or Boost
///        License as libdivide.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the zlib or Boost License. See the
/// COPYING file in the top level directory.
///

#ifndef FAST_DIV_128_64_HPP
#define FAST_DIV_128_64_HPP

#include <ctz.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <stdint.h>

#if defined(HAVE_INT128_T)

namespace primecount {

/// Prints an error message and calls std::abort()
[[noreturn]] void error_fast_div_128_64(uint128_t x, uint64_t y);

} // namespace

namespace {

using primecount::uint128_t;

/// primecount's performance depends heavily on the speed of the
/// 128-bit / 64-bit = 64-bit integer division. When using GCC's and
/// Clang's 128-bit integer division implementations this usually
/// emits a function call to __udivti3(a, b). According to my
/// benchmarks this function call adds significant overhead and hence
/// we can improve performance by using our own implementation
/// that is guaranteed to be inlined.
///
/// This algorithm is a modified version of libdivide's
/// libdivide_128_div_64_to_64() function which is heavily optimized
/// and very fast in practice.
///
ALWAYS_INLINE uint64_t fast_div_128_64(uint128_t x, uint64_t y)
{
  ASSERT(y > 0);
  uint64_t numlo = uint64_t(x);
  uint64_t numhi = uint64_t(x >> 64);
  uint64_t den = y;

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))
  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(numlo), "+d"(numhi) : [divider] "r"(den));

  return numlo;
#else
  // Use 64-bit integer division if possible.
  if (numhi == 0)
    return numlo / den;

  // Check for 64-bit overflow of quotient.
  if_unlikely(numhi >= den)
    primecount::error_fast_div_128_64(x, den);

  // We work in base 2**32.
  // A uint32 holds a single digit. A uint64 holds two digits.
  // Our numerator is conceptually [num3, num2, num1, num0].
  // Our denominator is [den1, den0].
  constexpr uint64_t b = uint64_t(1) << 32;

  // The high and low digits of our computed quotient.
  uint32_t q1, q0;
  // The normalization shift factor.
  int shift;
  // The high and low digits of our denominator (after normalizing).
  // Also the low 2 digits of our numerator (after normalizing).
  uint32_t den1, den0, num1, num0;
  // A partial remainder.
  uint64_t rem;
  // The estimated quotient, and its corresponding
  // remainder (unrelated to true remainder).
  uint64_t qhat, rhat;
  // Variables used to correct the estimated quotient.
  uint64_t c1, c2;

  // Determine the normalization factor. We multiply den by this, so that its leading digit is at
  // least half b. In binary this means just shifting left by the number of leading zeros,
  // so that there's a 1 in the MSB.
  // We also shift numer by the same amount. This cannot overflow because numhi < den.
  // The expression (-shift & 63) is the same as (64 - shift), except it avoids the UB of shifting
  // by 64. The funny bitwise 'and' ensures that numlo does not get shifted into numhi if shift is
  // 0. clang 11 has an x86 codegen bug here: see LLVM bug 50118. The sequence below avoids it.
  shift = clz64(den);
  den <<= shift;
  numhi <<= shift;
  numhi |= (numlo >> (-shift & 63)) & uint64_t(-int64_t(shift) >> 63);
  numlo <<= shift;

  // Extract the low digits of the numerator and
  // both digits of the denominator.
  num1 = uint32_t(numlo >> 32);
  num0 = uint32_t(numlo);
  den1 = uint32_t(den >> 32);
  den0 = uint32_t(den);

  // We wish to compute q1 = [n3 n2 n1] / [d1 d0].
  // Estimate q1 as [n3 n2] / [d1], and then correct it.
  // Note while qhat may be 2 digits, q1 is always 1 digit.
  qhat = numhi / den1;
  rhat = numhi % den1;
  c1 = qhat * den0;
  c2 = rhat * b + num1;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) + 1;
  q1 = uint32_t(qhat);

  // Compute the true (partial) remainder.
  rem = numhi * b + num1 - q1 * den;

  // We wish to compute q0 = [rem1 rem0 n0] / [d1 d0].
  // Estimate q0 as [rem1 rem0] / [d1] and correct it.
  qhat = rem / den1;
  rhat = rem % den1;
  c1 = qhat * den0;
  c2 = rhat * b + num0;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) + 1;
  q0 = uint32_t(qhat);

  uint64_t q = (uint64_t(q1) << 32) | q0;
  ASSERT(q == x / y);

  return q;
#endif
}

} // namespace

#endif
#endif
