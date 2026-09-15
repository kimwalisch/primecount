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

#include <cpu_arch_macros.hpp>
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

#if defined(HAVE_INT128_T)

/// Defined in src/udiv128.hpp
uint64_t udiv_128_by_64_to_64(uint64_t hi, uint64_t lo, uint64_t d);
uint128_t udiv_128_by_64_to_128(uint64_t hi, uint64_t lo, uint64_t d);
uint128_t udiv_128_by_128_to_128(uint64_t hi, uint64_t lo, uint64_t d1, uint64_t d0);

/// Used for (128-bit / 32-bit) = 64-bit.
/// This is Knuth's short division by a single-precision integer,
/// specialized to base 2^32 and a quotient that fits into 64 bits.
/// See The Art of Computer Programming, Volume 2, Section 4.3.1,
/// Exercise 16.
///
ALWAYS_INLINE uint64_t udiv_128_by_32_to_64(uint64_t hi,
                                            uint64_t lo,
                                            uint32_t d)
{
  ASSERT(d != 0);
  ASSERT(hi < d);

  // If the numerator fits in 64 bits,
  // use a native 64-bit division.
  if (hi == 0)
    return lo / d;

  uint64_t n1 = (hi << 32) | (lo >> 32);
  uint64_t n0 = lo;
  uint64_t dhi = uint64_t(d) << 32;
  uint64_t q1 = n1 / d;
  uint64_t dividend = n0 - q1 * dhi;
  uint64_t q0 = dividend / d;

  return (q1 << 32) | q0;
}

#endif

/// Used for (64-bit / 32-bit) = 64-bit
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint32_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

#if defined(ENABLE_DIV32) && \
    defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  uint32_t hi = uint32_t(uint64_t(x) >> 32);

  if (hi < uint32_t(y))
  {
    uint32_t lo = uint32_t(x);
    uint32_t d = y;

    // (64-bit / 32-bit) = 32-bit.
    // When we know the result fits into 32-bit (even
    // though the numerator is 64-bit) we can use the divl
    // instruction instead of doing a full 64-bit division.
    __asm__("divl %[divider]"
            : "+a"(lo), "+d"(hi) : [divider] "r"(d));

    return lo;
  }
#endif

  // Unsigned integer division is usually
  // faster than signed integer division.
  return uint64_t(x) / uint32_t(y);
}

/// Used for (64-bit / 64-bit) = 64-bit
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                                       sizeof(Y) == sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  return uint64_t(x) / uint64_t(y);
}

#if defined(HAVE_INT128_T)

/// Used for (128-bit / 32-bit) = 128-bit.
/// Used for (128-bit / 64-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint128_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  uint64_t hi = uint64_t(uint128_t(x) >> 64);
  uint64_t lo = uint64_t(x);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  if (hi < uint64_t(y))
  {
    uint64_t lo = uint64_t(x);
    uint64_t d = y;

    // (128-bit / 64-bit) = 64-bit.
    // When we know the result fits into 64-bit (even
    // though the numerator is 128-bit) we can use the divq
    // instruction instead of doing a full 128-bit division.
    __asm__("div %[divider]"
            : "+a"(lo), "+d"(hi) : [divider] "r"(d));

    return lo;
  }
#endif

  return udiv_128_by_64_to_128(hi, lo, y);
}

/// Used for (128-bit / 128-bit) = 128-bit
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint128_t) &&
                                       sizeof(Y) == sizeof(uint128_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  uint64_t hi = uint64_t(uint128_t(x) >> 64);
  uint64_t lo = uint64_t(x);
  uint64_t d1 = uint64_t(uint128_t(y) >> 64);
  uint64_t d0 = uint64_t(y);

  return udiv_128_by_128_to_128(hi, lo, d1, d0);
}

#endif

/// Used for (64-bit / 32-bit) = 64-bit.
/// Used for (64-bit / 64-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(X)), uint64_t>::type
fast_div64(X x, Y y)
{
  return (uint64_t) fast_div(x, y);
}

#if defined(HAVE_INT128_T)

/// Used for (128-bit / 32-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint128_t) &&
                                       sizeof(Y) <= sizeof(uint32_t)), uint64_t>::type
fast_div64(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  uint64_t hi = uint64_t(uint128_t(x) >> 64);
  uint64_t lo = uint64_t(x);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  uint64_t d = y;
  ASSERT(hi < d);

  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(lo), "+d"(hi) : [divider] "r"(d));

  return lo;
#else
  return udiv_128_by_32_to_64(hi, lo, y);
#endif
}

/// Used for (128-bit / 64-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint128_t) &&
                                       sizeof(Y) == sizeof(uint64_t)), uint64_t>::type
fast_div64(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  uint64_t hi = uint64_t(uint128_t(x) >> 64);
  uint64_t lo = uint64_t(x);
  uint64_t d = y;

  ASSERT(hi < d);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(lo), "+d"(hi) : [divider] "r"(d));

  return lo;
#else
  return udiv_128_by_64_to_64(hi, lo, d);
#endif
}

#endif

#if defined(HAVE_INT128_T) && \
   (defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE))

/// Used for (128-bit / 32-bit) = 64-bit.
/// This is Knuth's short division by a single-precision integer,
/// specialized to base 2^32 and a quotient that fits into 64 bits.
/// The remainder update is rearranged to shorten the SVE dependency chain.
/// See The Art of Computer Programming, Volume 2, Section 4.3.1,
/// Exercise 16.
///
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div_128_by_32_to_64(svbool_t pg,
                                                 uint128_t numer,
                                                 svuint64_t divisor)
{
  uint64_t num1 = uint64_t(numer >> 32);
  uint64_t num0 = uint64_t(numer);

  svuint64_t denhi = svlsl_n_u64_x(pg, divisor, 32);
  svuint64_t dividend = svdup_n_u64(num1);
  svuint64_t q1 = svdiv_u64_x(pg, dividend, divisor);
  dividend = svdup_n_u64(num0);
  dividend = svmls_u64_x(pg, dividend, q1, denhi);
  svuint64_t q0 = svdiv_u64_x(pg, dividend, divisor);

  return svorr_u64_x(pg, svlsl_n_u64_x(pg, q1, 32), q0);
}

/// Used for (128-bit / 64-bit) = 64-bit.
/// This is the branchless correction variant of Knuth Algorithm D used by
/// libdivide's divllu() implementation, vectorized across the divisors.
/// The numerator is common to all lanes while each lane has its own divisor.
/// https://github.com/ridiculousfish/libdivide/blob/master/doc/divlu.c
///
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div_128_by_64_to_64(svbool_t pg,
                                                 uint128_t numer,
                                                 svuint64_t divisor)
{
  uint64_t numer_lo = uint64_t(numer);
  uint64_t numer_hi = uint64_t(numer >> 64);

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
  svuint64_t den0 = svand_n_u64_x(pg, den, 0xffffffff);
  svuint64_t denhi = svand_n_u64_x(pg, den, 0xffffffff00000000);
  svuint64_t num1 = svlsr_n_u64_x(pg, numlo, 32);
  svuint64_t num0 = svand_n_u64_x(pg, numlo, 0xffffffff);

  // Estimate and correct the high 32 quotient bits.
  // Form c2 = (numhi - q1 * den1) * 2^32 + num1 modulo 2^64.
  // Shift and add before MLS so they do not depend on the quotient.
  svuint64_t q1 = svdiv_u64_x(pg, numhi, den1);
  svuint64_t c1 = svmul_u64_x(pg, q1, den0);
  svuint64_t c2 = svlsl_n_u64_x(pg, numhi, 32);
  c2 = svadd_u64_x(pg, c2, num1);
  c2 = svmls_u64_x(pg, c2, q1, denhi);

  svbool_t corr1 = svcmpgt_u64(pg, c1, c2);
  svuint64_t delta = svsub_u64_x(pg, c1, c2);
  svbool_t corr2 = svcmpgt_u64(corr1, delta, den);
  q1 = svsub_n_u64_m(corr1, q1, 1);
  q1 = svsub_n_u64_m(corr2, q1, 1);

  // True partial remainder needed to estimate the low 32 quotient bits.
  // rem = c2 - c1 + correction * den.
  // Reuse the correction predicates, avoiding q1 * den multiply.
  svuint64_t rem = svsub_u64_x(pg, c2, c1);
  rem = svadd_u64_m(corr1, rem, den);
  rem = svadd_u64_m(corr2, rem, den);

  // Estimate and correct the low 32 quotient bits.
  svuint64_t q0 = svdiv_u64_x(pg, rem, den1);
  c1 = svmul_u64_x(pg, q0, den0);
  c2 = svlsl_n_u64_x(pg, rem, 32);
  c2 = svadd_u64_x(pg, c2, num0);
  c2 = svmls_u64_x(pg, c2, q0, denhi);

  corr1 = svcmpgt_u64(pg, c1, c2);
  delta = svsub_u64_x(pg, c1, c2);
  corr2 = svcmpgt_u64(corr1, delta, den);
  q0 = svsub_n_u64_m(corr1, q0, 1);
  q0 = svsub_n_u64_m(corr2, q0, 1);

  return svorr_u64_x(pg, svlsl_n_u64_x(pg, q1, 32), q0);
}

// -------------------------------------------------------------------

/// Used for (64-bit / 64-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint64_t numer,
                                   svuint64_t divisor)
{
  return svdiv_u64_x(pg, svdup_n_u64(numer), divisor);
}

/// Used for (128-bit / 64-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint128_t numer,
                                   svuint64_t divisor)
{
  svbool_t all_divisors_64bit = svcmpgt_n_u64(pg, divisor, UINT32_MAX);

  // Use simpler base-2^32 long division
  // if all divisors fit into 32 bits.
  if (!svptest_any(pg, all_divisors_64bit))
    return sve_div_128_by_32_to_64(pg, numer, divisor);
  else
    return sve_div_128_by_64_to_64(pg, numer, divisor);
}

// -------------------------------------------------------------------

/// Used for (64-bit / 32-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint64_t numer,
                                   const uint32_t* divisors)
{
  svuint64_t divisor = svld1uw_u64(pg, divisors);
  return svdiv_u64_x(pg, svdup_n_u64(numer), divisor);
}

/// Used for (64-bit / 64-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint64_t numer,
                                   const int64_t* divisors)
{
  svuint64_t divisor = svreinterpret_u64_s64(svld1_s64(pg, divisors));
  return svdiv_u64_x(pg, svdup_n_u64(numer), divisor);
}

/// Used for (128-bit / 32-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint128_t numer,
                                   const uint32_t* divisors)
{
  svuint64_t divisor = svld1uw_u64(pg, divisors);
  return sve_div_128_by_32_to_64(pg, numer, divisor);
}

/// Used for (128-bit / 64-bit) = 64-bit
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
ALWAYS_INLINE svuint64_t sve_div64(svbool_t pg,
                                   uint128_t numer,
                                   const int64_t* divisors)
{
  svuint64_t divisor = svreinterpret_u64_s64(svld1_s64(pg, divisors));
  return sve_div64(pg, numer, divisor);
}

// -------------------------------------------------------------------

#endif

} // namespace

#endif
