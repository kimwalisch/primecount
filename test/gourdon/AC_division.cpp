///
/// @file  AC_division.cpp
/// @brief Test AC ARM SVE division, masked tails, and pi accumulation.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <cpu_arch_macros.hpp>
#include <fast_div.hpp>
#include <int128_t.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <LoadBalancerAC.hpp>
#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>

#if defined(ENABLE_ARM_SVE)
  #include "AC_arm_sve.hpp"
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #include "AC_arm_sve.hpp"
  #include <cpu_supports_arm_sve.hpp>
#endif

using namespace primecount;

namespace {

MAYBE_UNUSED void check(bool ok)
{
  if (!ok)
  {
    std::cerr << "AC ARM SVE division or masked pi lookup failed!" << std::endl;
    std::exit(1);
  }
}

#if defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE)

template <typename PRIME, typename XP>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_pi_arm_sve(XP xp)
{
  SegmentedPiTable segmentedPi;
  segmentedPi.init(128, 256, 256);

  for (uint64_t size = 1; size <= 65; size++)
  {
    Vector<PRIME> primes(size + 1);
    primes[0] = 0;
    for (uint64_t i = 1; i <= size; i++)
      primes[i] = (PRIME) (xp / (128 + i));

    for (uint64_t start = 1; start <= size + 1; start++)
    {
      uint64_t expected1 = 0;
      uint64_t expected2 = 0;
      for (uint64_t i = start; i <= size; i++)
      {
        uint64_t quotient = uint64_t(xp / primes[i]);
        uint64_t count = segmentedPi[quotient];
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
  #if defined(ENABLE_ARM_SVE)
    check_pi_arm_sve<uint32_t>(uint64_t(1) << 20);
    check_pi_arm_sve<int64_t>(uint64_t(1) << 40);
    check_pi_arm_sve<int64_t>(UINT64_MAX);

    #if defined(HAVE_INT128_T)
      check_pi_arm_sve<int64_t>(uint128_t(1) << 68);
    #endif

    std::cout << "ARM SVE division and masked pi lookups passed." << std::endl;

  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
    {
      check_pi_arm_sve<uint32_t>(uint64_t(1) << 20);
      check_pi_arm_sve<int64_t>(uint64_t(1) << 40);
      check_pi_arm_sve<int64_t>(UINT64_MAX);

      #if defined(HAVE_INT128_T)
        check_pi_arm_sve<int64_t>(uint128_t(1) << 68);
      #endif

      std::cout << "ARM SVE division and masked pi lookups passed." << std::endl;
    }
  #else
    std::cout << "SIMD execution skipped: no supported AC SIMD backend on this CPU/build." << std::endl;
  #endif

  return 0;
}
