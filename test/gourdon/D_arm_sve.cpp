///
/// @file  D_arm_sve.cpp
/// @brief Test ARM SVE factor index conversion and batch division.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <BaseFactorTable.hpp>
#include <cpu_arch_macros.hpp>
#include <fast_div.hpp>
#include <int128_t.hpp>
#include <LoadBalancerS2.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <phi_vector.hpp>
#include <sieve/Sieve.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>

#if defined(ENABLE_ARM_SVE)
  #include "D_arm_sve.hpp"
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #include "D_arm_sve.hpp"
  #include <cpu_supports_arm_sve.hpp>
#endif

using namespace primecount;

namespace {

MAYBE_UNUSED void check(bool ok)
{
  if (!ok)
  {
    std::cerr << "ARM SVE factor conversion or batch division failed!" << std::endl;
    std::exit(1);
  }
}

#if defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE)

template <typename Index>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_conversion_arm_sve(uint64_t base)
{
  Vector<Index> indexes(480);
  Vector<uint64_t> numbers(480);

  for (uint64_t i = 0; i < indexes.size(); i++)
    indexes[i] = Index(base + i);

  for (uint64_t size = 1; size <= indexes.size(); size++)
  {
    for (uint64_t i = 0; i < size; i += svcntd())
    {
      svbool_t pg = svwhilelt_b64(i, size);
      svuint64_t m = BaseFactorTable::to_number_arm_sve(pg, &indexes[i]);
      svst1_u64(pg, &numbers[i], m);
    }

    for (uint64_t i = 0; i < size; i++)
      check(numbers[i] == uint64_t(BaseFactorTable::to_number(indexes[i])));
  }
}

template <typename Index>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void check_batch_arm_sve(uint64_t base)
{
  Array<Index, 128> indexes;
  Array<int64_t, 129> results;

  for (uint64_t i = 0; i < indexes.size(); i++)
    indexes[i] = Index(base + indexes.size() - 1 - i);

  uint64_t min_m = BaseFactorTable::to_number(base);
  const uint64_t numerators[] = { 0, min_m - 1, min_m, UINT64_MAX };

  for (uint64_t size = 0; size <= indexes.size(); size++)
  {
    for (uint64_t xp : numerators)
    {
      results[size] = -1;
      batch_div_arm_sve(xp, indexes, results, size, 0);
      check(results[size] == -1);

      for (uint64_t i = 0; i < size; i++)
        check(results[i] == int64_t(xp / BaseFactorTable::to_number(indexes[i])));
    }

    // Test the low subtraction added to batch_div_arm_sve().
    {
      uint64_t xp = UINT64_MAX;
      uint64_t low = 1;
      results[size] = -1;
      batch_div_arm_sve(xp, indexes, results, size, low);
      check(results[size] == -1);

      for (uint64_t i = 0; i < size; i++)
      {
        uint64_t m = BaseFactorTable::to_number(indexes[i]);
        check(results[i] == int64_t(xp / m - low));
      }
    }

    #ifdef HAVE_INT128_T
      uint128_t xp = uint128_t(min_m) * INT64_MAX;

      for (uint64_t remainder : { uint64_t(0), min_m - 1 })
      {
        results[size] = -1;
        batch_div_arm_sve(xp + remainder, indexes, results, size, 0);
        check(results[size] == -1);

        for (uint64_t i = 0; i < size; i++)
          check(results[i] == int64_t((xp + remainder) / BaseFactorTable::to_number(indexes[i])));
      }
    #endif
  }
}

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
void test_arm_sve()
{
  uint64_t max_index32 = UINT32_MAX;
  uint64_t max_index64 = BaseFactorTable::to_index(INT64_MAX);

  check_conversion_arm_sve<uint32_t>(0);
  check_conversion_arm_sve<uint32_t>(max_index32 / 480 * 480 - 480);
  check_conversion_arm_sve<uint32_t>(max_index32 - 479);
  check_conversion_arm_sve<int64_t>(0);
  check_conversion_arm_sve<int64_t>(max_index32 - 239);
  check_conversion_arm_sve<int64_t>(max_index64 - 479);

  check_batch_arm_sve<uint32_t>(1);
  check_batch_arm_sve<uint32_t>(BaseFactorTable::to_index(1ull << 32) - 64);
  check_batch_arm_sve<uint32_t>(max_index32 - 127);
  check_batch_arm_sve<int64_t>(1);
  check_batch_arm_sve<int64_t>(BaseFactorTable::to_index(1ull << 32) - 64);
  check_batch_arm_sve<int64_t>(max_index32 - 63);
  check_batch_arm_sve<int64_t>(max_index64 - 127);

  std::cout << "ARM SVE factor conversion and batch division passed." << std::endl;
}

#endif

} // namespace

int main()
{
  #if defined(ENABLE_ARM_SVE)
    test_arm_sve();
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      test_arm_sve();
  #else
    std::cout << "ARM SVE tests skipped: no supported backend in this build." << std::endl;
  #endif

  return 0;
}
