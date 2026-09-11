///
/// @file  fast_div.hpp
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

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <macros.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <type_traits>

#if defined(HAVE_INT128_T) && \
   (defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE))
  #include <arm_sve.h>
#endif

namespace primecount {

/// Used for (64-bit / 32-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint32_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;

#if defined(ENABLE_DIV32) && \
    defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  uint32_t high = uint32_t(UX(x) >> 32);

  if (high < UY(y))
  {
    uint32_t low = uint32_t(x);
    uint32_t d = y;

    // (64-bit / 32-bit) = 32-bit.
    // When we know the result fits into 32-bit (even
    // though the numerator is 64-bit) we can use the divl
    // instruction instead of doing a full 64-bit division.
    __asm__("divl %[divider]"
            : "+a"(low), "+d"(high) : [divider] "r"(d));

    return low;
  }
#endif

  return UX(x) / UY(y);
}

/// Used for (128-bit / 32-bit) = 128-bit.
/// Used for (128-bit / 64-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  uint64_t high = uint64_t(UX(x) >> 64);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  if (high < UY(y))
  {
    uint64_t low = uint64_t(x);
    uint64_t d = y;

    // (128-bit / 64-bit) = 64-bit.
    // When we know the result fits into 64-bit (even
    // though the numerator is 128-bit) we can use the divq
    // instruction instead of doing a full 128-bit division.
    __asm__("div %[divider]"
            : "+a"(low), "+d"(high) : [divider] "r"(d));

    return low;
  }
#else
  // This optimization is very important on non x64 CPUs
  // such as ARM64. On AWS Graviton 4 CPUs it e.g. improves
  // performance by about 60% when computing AC(1e22).
  if (high == 0)
    return uint64_t(x) / UY(y);
#endif

  return UX(x) / UY(y);
}

/// Used for  (64-bit /  64-bit) =  64-bit.
/// Used for (128-bit / 128-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) >= sizeof(uint64_t) &&
                                       sizeof(Y) == sizeof(X)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  return UX(x) / UY(y);
}

/// Used for (128-bit / 32-bit) = 64-bit.
/// Used for (128-bit / 64-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), uint64_t>::type
fast_div64(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  using UX = typename pstd::make_unsigned<X>::type;

  uint64_t low = uint64_t(x);
  uint64_t high = uint64_t(UX(x) >> 64);
  uint64_t d = y;

  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(low), "+d"(high) : [divider] "r"(d));

  return low;
#else
  return (uint64_t) fast_div(x, y);
#endif
}

/// Used for (64-bit / 32-bit) = 64-bit.
/// Used for (64-bit / 64-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(X)), uint64_t>::type
fast_div64(X x, Y y)
{
  return (uint64_t) fast_div(x, y);
}

#if defined(HAVE_INT128_T) && \
   (defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE))

/// Narrowing unsigned division: uint128_t / uint64_t -> uint64_t.
/// The quotient must fit into 64 bits, i.e. (numer >> 64) < divisor.
///
/// This is the branchless correction variant of Knuth Algorithm D used by
/// libdivide's divllu() implementation, vectorized across the divisors.
/// The numerator is common to all lanes while each lane has its own divisor.
/// https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t fast_div64(svbool_t pg,
                                    uint128_t numer,
                                    svuint64_t divisor)
{
  const uint64_t numer_lo = uint64_t(numer);
  const uint64_t numer_hi = uint64_t(numer >> 64);

  svuint64_t shift = svclz_u64_x(pg, divisor);
  svuint64_t den = svlsl_u64_x(pg, divisor, shift);

  svuint64_t lo = svdup_n_u64(numer_lo);
  svuint64_t hi = svdup_n_u64(numer_hi);
  svuint64_t inv_shift = svsub_u64_x(pg, svdup_n_u64(64), shift);

  // Normalize the 128-bit numerator. SVE shifts use the full shift count,
  // hence shifting right by 64 when shift == 0 yields zero as required.
  svuint64_t numhi = svlsl_u64_x(pg, hi, shift);
  numhi = svorr_u64_x(pg, numhi, svlsr_u64_x(pg, lo, inv_shift));
  svuint64_t numlo = svlsl_u64_x(pg, lo, shift);

  svuint64_t den1 = svlsr_n_u64_x(pg, den, 32);
  svuint64_t den0 = svand_n_u64_x(pg, den, 0xffffffffu);
  svuint64_t num1 = svlsr_n_u64_x(pg, numlo, 32);
  svuint64_t num0 = svand_n_u64_x(pg, numlo, 0xffffffffu);

  // Estimate and correct the high 32 quotient bits.
  svuint64_t q1 = svdiv_u64_x(pg, numhi, den1);
  svuint64_t rhat = svmls_u64_x(pg, numhi, q1, den1);
  svuint64_t c1 = svmul_u64_x(pg, q1, den0);
  svuint64_t c2 = svlsl_n_u64_x(pg, rhat, 32);
  c2 = svadd_u64_x(pg, c2, num1);

  svbool_t corr1 = svcmpgt_u64(pg, c1, c2);
  svuint64_t delta = svsub_u64_x(pg, c1, c2);
  svbool_t corr2 = svcmpgt_u64(corr1, delta, den);
  q1 = svsub_n_u64_m(corr1, q1, 1);
  q1 = svsub_n_u64_m(corr2, q1, 1);

  // True partial remainder needed to estimate the low 32 quotient bits.
  // From numhi = qhat * den1 + rhat:
  // rem = c2 - c1 + correction * den.
  // Reuse the correction predicates, avoiding q1 * den multiply.
  svuint64_t rem = svsub_u64_x(pg, c2, c1);
  rem = svadd_u64_m(corr1, rem, den);
  rem = svadd_u64_m(corr2, rem, den);

  // Estimate and correct the low 32 quotient bits.
  svuint64_t q0 = svdiv_u64_x(pg, rem, den1);
  rhat = svmls_u64_x(pg, rem, q0, den1);
  c1 = svmul_u64_x(pg, q0, den0);
  c2 = svlsl_n_u64_x(pg, rhat, 32);
  c2 = svadd_u64_x(pg, c2, num0);

  corr1 = svcmpgt_u64(pg, c1, c2);
  delta = svsub_u64_x(pg, c1, c2);
  corr2 = svcmpgt_u64(corr1, delta, den);
  q0 = svsub_n_u64_m(corr1, q0, 1);
  q0 = svsub_n_u64_m(corr2, q0, 1);

  return svorr_u64_x(pg, svlsl_n_u64_x(pg, q1, 32), q0);
}

#endif

} // namespace

#endif
