///
/// @file   FactorTableD.cpp
/// @brief  FactorTableD is a compressed lookup table of mu (moebius),
///         lpf (least prime factor) and mpf (max prime factor).
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

// factorTable.mu(n) = 0 is disabled by default for performance
// reasons, we only enable it for testing.
#define ENABLE_MU_0_TESTING

#include <FactorTableD.hpp>
#include <generate_primes.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  check(FactorTableD<uint16_t>::max() == 4294705155ll);
  check(FactorTableD<uint32_t>::max() == pstd::numeric_limits<int64_t>::max());
  check(FactorTableDOrdinal<uint8_t>::max() == 2588880ll);
  check(FactorTableDOrdinal<uint16_t>::max() == 674982194328ll);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist_y(50000, 60000);
  std::uniform_int_distribution<int64_t> dist_z(1200000, 1500000);

  auto y = dist_y(gen);
  auto z = dist_z(gen);
  auto threads = get_num_threads();
  auto primes = generate_primes<uint32_t>(z);
  auto lpf = generate_lpf(z, primes);
  auto mpf = generate_mpf(z, primes);
  auto mu = generate_moebius(z, primes);
  Vector<uint32_t> prime_indexes(z + 1);

  for (std::size_t i = 1; i < primes.size(); i++)
    prime_indexes[primes[i]] = (uint32_t) i;

  FactorTableD<uint16_t> factorTable(y, z, threads);
  FactorTableDOrdinal<uint16_t> ordinalTable(y, z, threads);
  int64_t uint16_max = pstd::numeric_limits<uint16_t>::max();
  int64_t limit = factorTable.first_coprime();
  std::vector<int> filter_primes = { 13, 17, 101, 9973 };

  for (int64_t n = 1; n <= z; n++)
  {
    int64_t i = factorTable.to_index(n);
    bool is_prime = (lpf[n] == n);

    // Check if n is coprime to the primes < limit
    for (int64_t j = 1; primes[j] < limit; j++)
      if (n % primes[j] == 0)
        goto not_coprime;

    // primes > y and square free numbers with a prime factor > y
    // have been removed from the FactorTableD.
    if (mpf[n] > y)
    {
      std::cout << "prime_factor_larger_y(" << n << ") = " << (factorTable.is_leaf(i) == 0);
      check(factorTable.is_leaf(i) == 0);
      check(ordinalTable.is_leaf(i) == 0);
      continue;
    }

    std::cout << "mu(" << n << ") = " << factorTable.mu(i);
    check(mu[n] == factorTable.mu(i));
    check(mu[n] == ordinalTable.mu(i));

    std::cout << "lpf(" << n << ") = " << lpf[n];

    // is_leaf(n) is a combination of the mu(n) (Möbius function),
    // lpf(n) (least prime factor) and mpf(n) (max prime factor)
    // functions. is_leaf(n) returns (with n = to_number(index)):
    //
    // 1) INT_MAX - 1  if n = 1
    // 2) INT_MAX      if n is a prime
    // 3) 0            if n has a prime factor > y
    // 4) 0            if moebius(n) = 0
    // 5) lpf - 1      if moebius(n) = 1
    // 6) lpf          if moebius(n) = -1

    if (n == 1)
    {
      check(factorTable.is_leaf(i) == uint16_max - 1);
      check(ordinalTable.is_leaf(i) == uint16_max - 1);
    }
    else if (is_prime)
    {
      check(factorTable.is_leaf(i) == uint16_max);
      check(ordinalTable.is_leaf(i) == uint16_max);
    }
    else if (mu[n] == 0)
    {
      check(factorTable.is_leaf(i) == 0);
      check(ordinalTable.is_leaf(i) == 0);
    }
    else
    {
      check(lpf[n] == factorTable.is_leaf(i) + (factorTable.mu(i) == 1));
      check(prime_indexes[lpf[n]] == ordinalTable.is_leaf(i));
    }

    // The numerical and ordinal encodings must make the same
    // filtering decision for every prime threshold used by D.
    if (n != 1)
    {
      for (int prime : filter_primes)
      {
        int64_t prime_index = prime_indexes[prime];
        bool numerical = factorTable.is_leaf(i) > factorTable.get_filter_value(prime, prime_index);
        bool ordinal = ordinalTable.is_leaf(i) > ordinalTable.get_filter_value(prime, prime_index);
        check(numerical == ordinal);
      }
    }

    not_coprime:;
  }

  // Forced small-width boundary test. Codes 254 and 255 are
  // reserved for n = 1 and primes, hence ordinal 253 is the largest
  // least-prime ordinal that fits into uint8_t.
  {
    int64_t p253 = primesieve::nth_prime(253);
    int64_t p254 = primesieve::nth_prime(254);
    int64_t n = p253 * p254;
    FactorTableDOrdinal<uint8_t> smallOrdinalTable(FactorTableDOrdinal<uint8_t>::max(),
                                                   FactorTableDOrdinal<uint8_t>::max(),
                                                   threads);
    int64_t i = smallOrdinalTable.to_index(n);
    check(smallOrdinalTable.is_leaf(i) == 253);
    check(smallOrdinalTable.mu(i) == 1);
  }

  // Compare single-threaded and parallel construction. The larger z
  // crosses FactorTableD's threading threshold and exercises bitmap
  // word-aligned thread boundaries.
  {
    int64_t parallel_y = 1000000;
    int64_t parallel_z = 10000001;
    FactorTableDOrdinal<uint16_t> singleThread(parallel_y, parallel_z, 1);
    FactorTableDOrdinal<uint16_t> multiThread(parallel_y, parallel_z, 2);
    int64_t size = singleThread.to_index(parallel_z) + 1;
    bool equal = true;

    for (int64_t i = 0; equal && i < size; i++)
    {
      equal = singleThread.is_leaf(i) == multiThread.is_leaf(i);

      if (equal && singleThread.is_leaf(i) != 0)
        equal = singleThread.mu(i) == multiThread.mu(i);
    }

    check(equal);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
