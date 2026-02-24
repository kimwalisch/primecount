///
/// @file  fast_div_128_to_64.cpp
/// @brief Test fast_div_128_to_64(x, y) function
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <fast_div.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <limits>

#ifdef HAVE_INT128_T

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

#endif

int main()
{
#ifdef HAVE_INT128_T

  std::random_device rd;
  std::mt19937_64 gen(rd());

  // Denominator: choose from 1..UINT64_MAX (avoid zero)
  std::uniform_int_distribution<uint64_t> dist_u64(1, std::numeric_limits<uint64_t>::max());
  // Full-range 64-bit distribution for q and low parts
  std::uniform_int_distribution<uint64_t> dist_u64_full(0, std::numeric_limits<uint64_t>::max());

  {
    // Test numerator UINT64_MAX and denomiator 1
    uint128_t x = std::numeric_limits<uint64_t>::max();
    uint64_t res = fast_div_128_to_64(x, 1);
    std::cout << "fast_div_128_to_64(" << x << ", " << 1 << ") = " << res;
    check(res == x / 1);

    // Test denomiator UINT64_MAX
    uint64_t u64_max = std::numeric_limits<uint64_t>::max();
    x = uint128_t(u64_max) * u64_max - 1;
    res = fast_div_128_to_64(x, u64_max);
    std::cout << "fast_div_128_to_64(" << x << ", " << u64_max << ") = " << res;
    check(res == u64_max - 1);

    x = uint128_t(u64_max) * u64_max;
    res = fast_div_128_to_64(x, u64_max);
    std::cout << "fast_div_128_to_64(" << x << ", " << u64_max << ") = " << res;
    check(res == u64_max);

    x = uint128_t(u64_max) * u64_max + (u64_max - 1);
    res = fast_div_128_to_64(x, u64_max);
    std::cout << "fast_div_128_to_64(" << x << ", " << u64_max << ") = " << res;
    check(res == u64_max);
  }

  // 2^n denominator tests
  for (int n = 0; n < 64; n++)
  {
    uint64_t den = uint64_t(1) << n;

    for (int j = 0; j < 10; j++)
    {
      uint64_t q = dist_u64_full(gen);

      // remainder r in [0, den-1]
      std::uniform_int_distribution<uint64_t> dist_r(0, den - 1);
      uint64_t r = dist_r(gen);
      // Build numerator: x = q * den + r
      uint128_t x = uint128_t(q) * uint128_t(den) + uint128_t(r);

      uint64_t res = fast_div_128_to_64(x, den);
      std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
      check(res == x / den);
    }
  }

  // 2^n+1 denominator tests
  for (int n = 1; n < 63; n++)
  {
    uint64_t den = (uint64_t(1) << n) + 1;

    for (int j = 0; j < 10; j++)
    {
      uint64_t q = dist_u64_full(gen);

      // remainder r in [0, den-1]
      std::uniform_int_distribution<uint64_t> dist_r(0, den - 1);
      uint64_t r = dist_r(gen);
      // Build numerator: x = q * den + r
      uint128_t x = uint128_t(q) * uint128_t(den) + uint128_t(r);

      uint64_t res = fast_div_128_to_64(x, den);
      std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
      check(res == x / den);
    }
  }

  // 2^n-1 denominator tests
  for (int n = 1; n < 64; n++)
  {
    uint64_t den = (uint64_t(1) << n) - 1;

    for (int j = 0; j < 10; j++)
    { 
      uint64_t q = dist_u64_full(gen);

      // remainder r in [0, den-1]
      std::uniform_int_distribution<uint64_t> dist_r(0, den - 1);
      uint64_t r = dist_r(gen);
      // Build numerator: x = q * den + r
      uint128_t x = uint128_t(q) * uint128_t(den) + uint128_t(r);

      uint64_t res = fast_div_128_to_64(x, den);
      std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
      check(res == x / den);
    }
  }

  // Test small numerators < 5000
  for (uint64_t x = 0; x < 5000; x++)
  {
    uint64_t den = dist_u64(gen) % 100 + 1;
    uint64_t res = fast_div_128_to_64(x, den);

    std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
    check(res == x / den);
  }

  const uint64_t edge_den[] = {
    1, 2, 3, 4, 5, std::numeric_limits<uint64_t>::max()
  };

  // Test denominator edge cases
  for (int i = 0; i < 1000; i++)
  {
    for (uint64_t den : edge_den)
    {
      uint64_t q = dist_u64_full(gen);

      // remainder r in [0, den-1]
      std::uniform_int_distribution<uint64_t> dist_r(0, den - 1);
      uint64_t r = dist_r(gen);
      // Build numerator: x = q * den + r
      uint128_t x = uint128_t(q) * uint128_t(den) + uint128_t(r);

      uint64_t res = fast_div_128_to_64(x, den);
      std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
      check(res == q);
    }
  }

  // Test random 64-bit numerator and denominator
  for (int i = 0; i < 3000; i++)
  {
    uint64_t x = dist_u64_full(gen);
    uint64_t den = dist_u64(gen);
    uint64_t res = fast_div_128_to_64(x, den);

    std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
    check(res == x / den);
  }

  // Test random 128-bit numerator and 64-bit denominator
  for (int i = 0; i < 7000; i++)
  {
    uint64_t den = dist_u64(gen);
    uint64_t q = dist_u64_full(gen);

    // remainder r in [0, den-1]
    std::uniform_int_distribution<uint64_t> dist_r(0, den - 1);
    uint64_t r = dist_r(gen);
    // Build numerator: x = q * den + r
    uint128_t x = uint128_t(q) * uint128_t(den) + uint128_t(r);

    uint64_t res = fast_div_128_to_64(x, den);
    std::cout << "fast_div_128_to_64(" << x << ", " << den << ") = " << res;
    check(res == q);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
