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

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("+sve")))
#endif
bool test_fast_div64_arm_sve(std::mt19937& gen)
{
  // ARM SVE supports at most 2048-bit vectors = 32 uint64_t lanes.
  constexpr uint64_t max_lanes = 32;
  uint64_t divisors[max_lanes] = {};
  uint64_t results[max_lanes] = {};
  uint64_t lanes = svcntd();

  if (lanes > max_lanes)
    return false;

  svbool_t all = svptrue_b64();

  // Test quotient boundary cases. The numerator is chosen as
  // divisor * UINT64_MAX + divisor - 1, hence the quotient is
  // exactly UINT64_MAX while still satisfying numer_hi < divisor.
  const uint64_t edge_divisors[] = {
    1, 2,
    uint64_t(1) << 31,
    uint64_t(1) << 32,
    (uint64_t(1) << 32) + 1,
    uint64_t(1) << 63,
    pstd::numeric_limits<uint64_t>::max()
  };

  for (uint64_t divisor : edge_divisors)
  {
    uint128_t numer = uint128_t(divisor) * pstd::numeric_limits<uint64_t>::max();
    numer += divisor - 1;

    svuint64_t den = svdup_n_u64(divisor);
    svuint64_t quot = fast_div64(all, numer, den);
    svst1_u64(all, results, quot);

    uint64_t expected = uint64_t(numer / divisor);
    for (uint64_t j = 0; j < lanes; j++)
      if (results[j] != expected)
        return false;
  }

  // Cases that exercise both 1-step and 2-step quotient corrections
  // in the underlying libdivide Algorithm-D formulation.
  struct TestCase
  {
    uint64_t hi;
    uint64_t lo;
    uint64_t divisor;
  };

  const TestCase correction_cases[] = {
    // q1: one correction.
    { 1663239473288121450ull, 8491104977468830630ull, 2789495100195680658ull },
    // q1 and q0: one correction each.
    { 5901740277908440862ull, 7753795759061867340ull, 15764240356045080578ull },
    // q1 and q0: two corrections each.
    { 5351022020821579747ull, 10990235090364358412ull, 9277324617043713967ull }
  };

  for (const TestCase& test : correction_cases)
  {
    uint128_t numer = (uint128_t(test.hi) << 64) | test.lo;
    svuint64_t den = svdup_n_u64(test.divisor);
    svuint64_t quot = fast_div64(all, numer, den);
    svst1_u64(all, results, quot);

    uint64_t expected = uint64_t(numer / test.divisor);
    for (uint64_t j = 0; j < lanes; j++)
      if (results[j] != expected)
        return false;
  }

  // Test random divisors with many different normalization shifts.
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
    uint128_t numer = (uint128_t(high) << 64) | low;

    for (uint64_t j = 0; j < lanes; j++)
    {
      // For wide numerators use bit widths 17..64. For 64-bit
      // numerators use 1..64, thereby exercising clz values 0..63.
      uint64_t bits = wide
        ? 17 + ((i + j) % 48)
        :  1 + ((i + j) % 64);
      uint64_t top_bit = uint64_t(1) << (bits - 1);
      divisors[j] = top_bit | (dist_u64(gen) & (top_bit - 1));
    }

    svuint64_t den = svld1_u64(all, divisors);
    svuint64_t quot = fast_div64(all, numer, den);
    svst1_u64(all, results, quot);

    for (uint64_t j = 0; j < lanes; j++)
      if (results[j] != uint64_t(numer / divisors[j]))
        return false;
  }

  // Test predication as used by the tail loop in AC_arm_sve.hpp.
  uint64_t active = lanes - 1;
  svbool_t pg = svwhilelt_b64(uint64_t(0), active);
  uint128_t numer = (uint128_t(UINT64_C(12345)) << 64) | UINT64_C(987654321);

  for (uint64_t j = 0; j < lanes; j++)
    divisors[j] = UINT64_C(65537) + j * 2;

  svuint64_t den = svld1_u64(pg, divisors);
  svuint64_t quot = fast_div64(pg, numer, den);
  svst1_u64(pg, results, quot);

  for (uint64_t j = 0; j < active; j++)
    if (results[j] != uint64_t(numer / divisors[j]))
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

  std::cout << "fast_div64(ARM SVE)";
  check(test_fast_div64_arm_sve(gen));

#elif defined(HAVE_INT128_T) && \
      defined(ENABLE_MULTIARCH_ARM_SVE)

  if (cpu_supports_sve)
  {
    std::cout << "fast_div64(ARM SVE)";
    check(test_fast_div64_arm_sve(gen));
  }

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
