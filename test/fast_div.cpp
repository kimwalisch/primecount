///
/// @file  fast_div.cpp
/// @brief Test fast_div(x, y) function
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <cpu_arch_macros.hpp>
#include <fast_div.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <random>

#if defined(HAVE_INT128_T) && \
    defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
#endif

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

#if defined(HAVE_INT128_T) && \
   (defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE))

constexpr uint64_t max_sve_lanes = 32;

template <typename Divisor>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
bool check_sve_div64(svbool_t pg,
                     uint64_t numer,
                     const Divisor* divisors)
{
  uint64_t results[max_sve_lanes];
  svuint64_t quot = sve_div64(pg, numer, divisors);
  svst1_u64(pg, results, quot);

  uint64_t active = svcntp_b64(pg, pg);
  for (uint64_t j = 0; j < active; j++)
    if (results[j] != numer / uint64_t(divisors[j]))
      return false;

  return true;
}

template <typename Divisor>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
bool check_sve_div64(svbool_t pg,
                     uint128_t numer,
                     const Divisor* divisors)
{
  uint64_t results[max_sve_lanes];
  svuint64_t quot = sve_div64(pg, numer, divisors);
  svst1_u64(pg, results, quot);

  uint64_t active = svcntp_b64(pg, pg);
  for (uint64_t j = 0; j < active; j++)
    if (results[j] != uint64_t(numer / uint64_t(divisors[j])))
      return false;

  return true;
}

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
bool test_sve_div64_arm_sve(std::mt19937& gen)
{
  // ARM SVE supports at most 2048-bit vectors = 32 uint64_t lanes.
  uint32_t divisors32[max_sve_lanes];
  int64_t divisors64[max_sve_lanes];
  uint64_t lanes = svcntd();

  if (lanes > max_sve_lanes)
    return false;

  svbool_t all = svptrue_b64();

  // Test quotient boundary cases for uint32_t divisors. The numerator
  // is chosen such that the quotient is exactly UINT64_MAX while still
  // satisfying numer_hi < divisor.
  const uint32_t edge_divisors32[] = {
    1, 2,
    uint32_t(1) << 31,
    pstd::numeric_limits<uint32_t>::max()
  };

  for (uint32_t divisor : edge_divisors32)
  {
    for (uint64_t j = 0; j < lanes; j++)
      divisors32[j] = divisor;

    uint128_t numer = uint128_t(divisor) * pstd::numeric_limits<uint64_t>::max();
    numer += divisor - 1;

    if (!check_sve_div64(all, numer, divisors32))
      return false;
  }

  // Test quotient boundary cases for int64_t divisors, including values
  // on both sides of the 2^32 boundary.
  const int64_t edge_divisors64[] = {
    1, 2,
    int64_t(1) << 31,
    pstd::numeric_limits<uint32_t>::max(),
    int64_t(1) << 32,
    (int64_t(1) << 32) + 1,
    int64_t(1) << 62,
    pstd::numeric_limits<int64_t>::max()
  };

  for (int64_t divisor : edge_divisors64)
  {
    for (uint64_t j = 0; j < lanes; j++)
      divisors64[j] = divisor;

    uint128_t numer = uint128_t(divisor) * pstd::numeric_limits<uint64_t>::max();
    numer += divisor - 1;

    if (!check_sve_div64(all, numer, divisors64))
      return false;
  }

  // Cases that exercise both 1-step and 2-step quotient corrections
  // in the underlying libdivide Algorithm-D formulation.
  struct TestCase
  {
    uint64_t hi;
    uint64_t lo;
    int64_t divisor;
  };

  const TestCase correction_cases[] = {
    // q1: one correction.
    { 1663239473288121450, 8491104977468830630, 2789495100195680658 },
    // q1 and q0: one correction each.
    { 9223372036854775806, UINT64_MAX, 9223372036854775807 },
    // q1 and q0: two corrections each.
    { 4611686022722355198, UINT64_MAX, 4611686022722355199 }
  };

  for (const TestCase& test : correction_cases)
  {
    for (uint64_t j = 0; j < lanes; j++)
      divisors64[j] = test.divisor;

    uint128_t numer = (uint128_t(test.hi) << 64) | test.lo;

    if (!check_sve_div64(all, numer, divisors64))
      return false;
  }

  // Test all four sve_div64() overload combinations using random divisors.
  // Every other iteration uses a true 128-bit numerator. In these
  // iterations high < 2^15 and divisor >= 2^16, which guarantees
  // that the quotient fits into uint64_t.
  std::uniform_int_distribution<uint64_t> dist_u64(
      0, pstd::numeric_limits<uint64_t>::max());

  for (uint64_t i = 0; i < 128; i++)
  {
    bool wide = i & 1;
    uint64_t high = wide ? 1 + (dist_u64(gen) & 0x7fff) : 0;
    uint64_t low = dist_u64(gen);
    uint128_t numer128 = (uint128_t(high) << 64) | low;
    uint64_t numer64 = dist_u64(gen);

    for (uint64_t j = 0; j < lanes; j++)
    {
      uint64_t bits32 = wide
        ? 17 + ((i + j) % 16)
        :  1 + ((i + j) % 32);
      uint64_t top_bit32 = uint64_t(1) << (bits32 - 1);
      divisors32[j] = uint32_t(top_bit32 |
          (dist_u64(gen) & (top_bit32 - 1)));

      uint64_t bits64 = wide
        ? 17 + ((i + j) % 47)
        :  1 + ((i + j) % 63);
      uint64_t top_bit64 = uint64_t(1) << (bits64 - 1);
      divisors64[j] = int64_t(top_bit64 |
          (dist_u64(gen) & (top_bit64 - 1)));
    }

    if (!check_sve_div64(all, numer64, divisors32) ||
        !check_sve_div64(all, numer64, divisors64) ||
        !check_sve_div64(all, numer128, divisors32) ||
        !check_sve_div64(all, numer128, divisors64))
      return false;
  }

  // Cover all normalization shifts, mixed divisor widths and large
  // numerators, including quotients near UINT64_MAX.
  for (uint64_t i = 0; i < 4096; i++)
  {
    uint64_t min_divisor = INT64_MAX;
    for (uint64_t j = 0; j < lanes; j++)
    {
      uint64_t top = 1ull << ((i + j) % 63);
      uint64_t divisor = top | (dist_u64(gen) & (top - 1));
      if (j == 0)
        divisor |= 1ull << 62;
      divisors64[j] = divisor;
      if (divisor < min_divisor)
        min_divisor = divisor;
    }

    uint64_t high = dist_u64(gen) % min_divisor;
    if (i & 1)
      high = min_divisor - 1;

    uint128_t numer = (uint128_t(high) << 64) | dist_u64(gen);
    uint64_t active = i % (lanes + 1);
    svbool_t pg = svwhilelt_b64(uint64_t(0), active);

    if (!check_sve_div64(all, numer, divisors64) ||
        !check_sve_div64(pg, numer, divisors64))
      return false;
  }

  // Test predication as used by the tail loop in AC_arm_sve.hpp.
  uint64_t active = lanes - 1;
  svbool_t pg = svwhilelt_b64(uint64_t(0), active);
  uint128_t numer128 = (uint128_t(12345) << 64) | 987654321;
  uint64_t numer64 = pstd::numeric_limits<uint64_t>::max();

  for (uint64_t j = 0; j < lanes; j++)
  {
    divisors32[j] = 65537 + uint32_t(j * 2);
    divisors64[j] = (int64_t(1) << 32) + 1 + int64_t(j * 2);
  }

  if (!check_sve_div64(pg, numer64, divisors32) ||
      !check_sve_div64(pg, numer64, divisors64) ||
      !check_sve_div64(pg, numer128, divisors32) ||
      !check_sve_div64(pg, numer128, divisors64))
    return false;

  return true;
}

#endif

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<int32_t> dist_i32(1, pstd::numeric_limits<int32_t>::max());
  std::uniform_int_distribution<uint64_t> dist_u64(0, pstd::numeric_limits<uint64_t>::max());

  // Test unsigned/signed
  for (int i = 0; i < 10000; i++)
  {
    uint64_t x = dist_i32(gen);
     int32_t y = dist_i32(gen);
    uint64_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_u64(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#ifdef HAVE_INT128_T

  std::uniform_int_distribution<uint64_t> dist_u62(0, uint64_t(1ull << 62));

  // Test signed/signed
  for (int i = 0; i < 10000; i++)
  {
    // Test x < 2^64
    int128_t x = dist_u64(gen);
     int32_t y = dist_i32(gen);
    int128_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    // Test x > 2^64
    int128_t low = dist_u64(gen);
    int128_t high = int128_t(dist_u62(gen)) << 64;
    x = high | low;
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#endif

#if defined(HAVE_INT128_T) && \
    defined(ENABLE_ARM_SVE)

  std::cout << "sve_div64(ARM SVE)";
  check(test_sve_div64_arm_sve(gen));

#elif defined(HAVE_INT128_T) && \
      defined(ENABLE_MULTIARCH_ARM_SVE)

  if (cpu_supports_sve)
  {
    std::cout << "sve_div64(ARM SVE)";
    check(test_sve_div64_arm_sve(gen));
  }

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
