///
/// @file  AC_division.cpp
/// @brief Test AC SIMD division, masked tails, and pi accumulation.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerAC.hpp>
#include <SegmentedPiTable.hpp>
#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <cpu_arch_macros.hpp>
#include <fast_div.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>

#if defined(ENABLE_LIBDIVIDE)
  #include "AC_libdivide.hpp"
  #if defined(ENABLE_AVX512_VPOPCNT) || \
      defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    #include "AC_libdivide_avx512.hpp"
    #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
      #include <cpu_supports_avx512_vpopcnt.hpp>
    #endif
  #elif defined(ENABLE_ARM_SVE) || \
        defined(ENABLE_MULTIARCH_ARM_SVE)
    #include "AC_libdivide_arm_sve.hpp"
    #if defined(ENABLE_MULTIARCH_ARM_SVE)
      #include <cpu_supports_arm_sve.hpp>
    #endif
  #endif
#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)
  #include "AC_arm_sve.hpp"
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    #include <cpu_supports_arm_sve.hpp>
  #endif
#endif

using namespace primecount;

namespace {

MAYBE_UNUSED void check(bool ok)
{
  if (!ok)
  {
    std::cerr << "AC vector division or masked pi lookup failed!" << std::endl;
    std::exit(1);
  }
}

#if defined(ENABLE_LIBDIVIDE)

#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
#endif
void check_vector(uint64_t xp,
                  const Vector<uint64_t>& primes,
                  const LibdividePrimes& lprimes)
{
  for (uint64_t i = 1; i < primes.size(); i += 8)
  {
    uint64_t count = min(uint64_t(8), primes.size() - i);
    __mmask8 mask = (__mmask8) (0xff >> (8 - count));
    Array<uint64_t, 8> quotients;
    quotients.fill(UINT64_MAX);
    __m512i numer = _mm512_set1_epi64(xp);
    __m512i q = divide_libdivide_avx512(numer, mask, &lprimes.magic[i], &lprimes.shift[i]);
    _mm512_mask_storeu_epi64(quotients.data(), mask, q);

    for (uint64_t j = 0; j < count; j++)
      check(quotients[j] == xp / primes[i + j]);

    for (uint64_t j = count; j < quotients.size(); j++)
      check(quotients[j] == UINT64_MAX);
  }
}

#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_vector(uint64_t xp,
                  const Vector<uint64_t>& primes,
                  const LibdividePrimes& lprimes)
{
  uint64_t lanes = svcntd();
  for (uint64_t i = 1; i < primes.size(); i += lanes)
  {
    uint64_t count = min(lanes, primes.size() - i);
    svbool_t pg = svwhilelt_b64(i, primes.size());
    Array<uint64_t, 32> quotients;
    quotients.fill(UINT64_MAX);
    svuint64_t numer = svdup_n_u64(xp);
    svuint64_t q = divide_libdivide_arm_sve(pg, numer, &lprimes.magic[i], &lprimes.shift[i]);
    svst1_u64(pg, quotients.data(), q);

    for (uint64_t j = 0; j < count; j++)
      check(quotients[j] == xp / primes[i + j]);

    for (uint64_t j = count; j < quotients.size(); j++)
      check(quotients[j] == UINT64_MAX);
  }
}

#endif

void check_division(uint64_t xp,
                    const Vector<uint64_t>& primes,
                    const LibdividePrimes& lprimes)
{
  for (uint64_t i = 1; i < primes.size(); i++)
    check(lprimes.divide(xp, i) == xp / primes[i]);

  #if defined(ENABLE_AVX512_VPOPCNT)
    check_vector(xp, primes, lprimes);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      check_vector(xp, primes, lprimes);
  #elif defined(ENABLE_ARM_SVE)
    check_vector(xp, primes, lprimes);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      check_vector(xp, primes, lprimes);
  #endif
}

void check_dividers()
{
  uint64_t random = 123456789;

  for (uint64_t size = 1; size <= 65; size++)
  {
    Vector<uint64_t> primes(size + 1);
    primes[0] = 0;

    for (uint64_t i = 1; i <= size; i++)
    {
      random ^= random << 13;
      random ^= random >> 7;
      random ^= random << 17;
      primes[i] = max(random, uint64_t(2));
    }

    // Include powers of two and divisors near the unsigned limit.
    primes[1] = 2;
    if (size > 1)
      primes[2] = UINT64_MAX;
    if (size > 2)
      primes[3] = uint64_t(1) << 63;
    if (size > 3)
      primes[4] = 3;

    LibdividePrimes lprimes(primes, 2);
    check_division(0, primes, lprimes);
    check_division(UINT64_MAX, primes, lprimes);
    check_division(random, primes, lprimes);

    for (int bit = 1; bit < 64; bit++)
    {
      uint64_t xp = uint64_t(1) << bit;
      check_division(xp - 1, primes, lprimes);
      check_division(xp, primes, lprimes);
      check_division(xp + 1, primes, lprimes);
    }
  }
}

#endif

#if defined(ENABLE_LIBDIVIDE) && \
    (defined(ENABLE_AVX512_VPOPCNT) || \
     defined(ENABLE_MULTIARCH_AVX512_VPOPCNT))

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
#endif
void check_pi_libdivide_avx512(uint64_t xp)
{
  SegmentedPiTable segmentedPi;
  segmentedPi.init(128, 256, 256);

  for (uint64_t size = 1; size <= 65; size++)
  {
    Vector<uint64_t> primes(size + 1);
    primes[0] = 0;
    for (uint64_t i = 1; i <= size; i++)
      primes[i] = xp / (128 + i);

    LibdividePrimes lprimes(primes, 1);

    for (uint64_t start = 1; start <= size + 1; start++)
    {
      uint64_t expected1 = 0;
      uint64_t expected2 = 0;
      for (uint64_t i = start; i <= size; i++)
      {
        uint64_t count = segmentedPi[xp / primes[i]];
        expected1 += count;
        expected2 += count * 2 - 7 + 2;
      }

      check(sum_pi_libdivide_avx512<uint64_t, 1>(xp, start, size, 2, lprimes, segmentedPi) == expected1);
      check(sum_pi_libdivide_avx512<uint64_t, 2>(xp, start, size, 7, lprimes, segmentedPi) == expected2);
      #ifdef HAVE_INT128_T
        check(sum_pi_libdivide_avx512<uint128_t, 1>(xp, start, size, 2, lprimes, segmentedPi) == expected1);
        check(sum_pi_libdivide_avx512<uint128_t, 2>(xp, start, size, 7, lprimes, segmentedPi) == expected2);
      #endif
    }
  }
}

#endif

#if defined(ENABLE_LIBDIVIDE) && \
    (defined(ENABLE_ARM_SVE) || \
     defined(ENABLE_MULTIARCH_ARM_SVE))

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_pi_libdivide_arm_sve(uint64_t xp)
{
  SegmentedPiTable segmentedPi;
  segmentedPi.init(128, 256, 256);

  for (uint64_t size = 1; size <= 65; size++)
  {
    Vector<uint64_t> primes(size + 1);
    primes[0] = 0;
    for (uint64_t i = 1; i <= size; i++)
      primes[i] = xp / (128 + i);

    LibdividePrimes lprimes(primes, 1);

    for (uint64_t start = 1; start <= size + 1; start++)
    {
      uint64_t expected1 = 0;
      uint64_t expected2 = 0;
      for (uint64_t i = start; i <= size; i++)
      {
        uint64_t count = segmentedPi[xp / primes[i]];
        expected1 += count;
        expected2 += count * 2 - 7 + 2;
      }

      check(sum_pi_libdivide_arm_sve<uint64_t, 1>(xp, start, size, 2, lprimes, segmentedPi) == expected1);
      check(sum_pi_libdivide_arm_sve<uint64_t, 2>(xp, start, size, 7, lprimes, segmentedPi) == expected2);
      #ifdef HAVE_INT128_T
        check(sum_pi_libdivide_arm_sve<uint128_t, 1>(xp, start, size, 2, lprimes, segmentedPi) == expected1);
        check(sum_pi_libdivide_arm_sve<uint128_t, 2>(xp, start, size, 7, lprimes, segmentedPi) == expected2);
      #endif
    }
  }
}

#endif

#if !defined(ENABLE_LIBDIVIDE) && \
    (defined(ENABLE_ARM_SVE) || \
     defined(ENABLE_MULTIARCH_ARM_SVE))

template <typename Prime>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_pi_arm_sve(uint64_t xp)
{
  SegmentedPiTable segmentedPi;
  segmentedPi.init(128, 256, 256);

  for (uint64_t size = 1; size <= 65; size++)
  {
    Vector<Prime> primes(size + 1);
    primes[0] = 0;
    for (uint64_t i = 1; i <= size; i++)
      primes[i] = xp / (128 + i);

    for (uint64_t start = 1; start <= size + 1; start++)
    {
      uint64_t expected1 = 0;
      uint64_t expected2 = 0;
      for (uint64_t i = start; i <= size; i++)
      {
        uint64_t count = segmentedPi[xp / primes[i]];
        expected1 += count;
        expected2 += count * 2 - 7 + 2;
      }

      check(sum_pi_arm_sve<uint64_t, 1>(xp, start, size, 2, primes, segmentedPi) == expected1);
      check(sum_pi_arm_sve<uint64_t, 2>(xp, start, size, 7, primes, segmentedPi) == expected2);
      #ifdef HAVE_INT128_T
        check(sum_pi_arm_sve<uint128_t, 1>(xp, start, size, 2, primes, segmentedPi) == expected1);
        check(sum_pi_arm_sve<uint128_t, 2>(xp, start, size, 7, primes, segmentedPi) == expected2);
      #endif
    }
  }
}

#endif

} // namespace

int main()
{
  #if defined(ENABLE_LIBDIVIDE)
    check_dividers();
    std::cout << "Decomposed libdivide parameters passed." << std::endl;
  #endif

  #if defined(ENABLE_LIBDIVIDE) && \
      (defined(ENABLE_AVX512_VPOPCNT) || \
       defined(ENABLE_MULTIARCH_AVX512_VPOPCNT))
    #if !defined(ENABLE_AVX512_VPOPCNT)
      if (cpu_supports_avx512_vpopcnt)
    #endif
    {
      check_pi_libdivide_avx512(uint64_t(1) << 40);
      std::cout << "AVX512 division and masked pi lookups passed." << std::endl;
      return 0;
    }
  #elif defined(ENABLE_ARM_SVE) || \
        defined(ENABLE_MULTIARCH_ARM_SVE)
    #if !defined(ENABLE_ARM_SVE)
      if (cpu_supports_sve)
    #endif
    {
      #if defined(ENABLE_LIBDIVIDE)
        check_pi_libdivide_arm_sve(uint64_t(1) << 40);
      #else
        check_pi_arm_sve<uint32_t>(uint64_t(1) << 20);
        check_pi_arm_sve<int64_t>(uint64_t(1) << 40);
      #endif
      std::cout << "ARM SVE division and masked pi lookups passed." << std::endl;
      return 0;
    }
  #endif

  std::cout << "SIMD execution skipped: no supported AC SIMD backend on this CPU/build." << std::endl;
  return 0;
}
