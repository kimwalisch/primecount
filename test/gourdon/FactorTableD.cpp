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
  if (!OK)
  {
    std::cout << "   ERROR\n";
    std::exit(1);
  }
}

int main()
{
  check(FactorTableD<uint16_t>::max() == 149060082888ll);
  check(FactorTableD<uint32_t>::max() == pstd::numeric_limits<int64_t>::max());
  check(FactorTableD<uint16_t>::to_index(
          FactorTableD<uint16_t>::max()) > UINT32_MAX);

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
  int64_t uint16_max = pstd::numeric_limits<uint16_t>::max();
  int64_t limit = factorTable.first_coprime();
  std::vector<uint32_t> filter_primes = { 13, 17, 101, 9973 };

  for (int64_t n = 1; n <= z; n++)
  {
    int64_t i = factorTable.to_index(n);
    bool is_prime = (lpf[n] == n);
    bool mu_OK;
    bool lpf_OK;

    // Check if n is coprime to the primes < limit
    for (int64_t j = 1; primes[j] < limit; j++)
      if (n % primes[j] == 0)
        goto not_coprime;

    // primes > y and square free numbers with a prime factor > y
    // have been removed from the FactorTableD.
    if (mpf[n] > y)
    {
      bool OK = (factorTable[i] == 0);

      if (!OK)
        std::cout << "prime_factor_larger_y(" << n << ") = " << factorTable[i];

      check(OK);
      continue;
    }

    mu_OK = (mu[n] == factorTable.mu(i));

    if (!mu_OK)
      std::cout << "mu(" << n << ") = " << factorTable.mu(i);

    check(mu_OK);

    // factorTable[index] is a combination of the mu(n) (Möbius
    // function), lpf(n) (least prime factor) and mpf(n) (max prime
    // factor) functions. It returns (with n = to_number(index)):
    //
    // 1) INT_MAX - 1    if n = 1
    // 2) INT_MAX        if n is a prime
    // 3) 0              if n has a prime factor > y
    // 4) 0              if moebius(n) = 0
    // 5) 2*pi(lpf)      if moebius(n) = 1
    // 6) 2*pi(lpf) + 1  if moebius(n) = -1

    if (n == 1)
      lpf_OK = (factorTable[i] == uint16_max - 1);
    else if (is_prime)
      lpf_OK = (factorTable[i] == uint16_max);
    else if (mu[n] == 0)
      lpf_OK = (factorTable[i] == 0);
    else
    {
      int64_t expected = prime_indexes[lpf[n]] * 2 + (mu[n] == -1);
      lpf_OK = (factorTable[i] == expected);
    }

    if (!lpf_OK)
      std::cout << "lpf(" << n << ") = " << lpf[n];

    check(lpf_OK);

    // Check that "factor[n] > encode(PrimePi(prime))" is equivalent
    // to "mu[n] != 0 && lpf[n] > prime" for a few prime thresholds.
    // (Numbers with mpf[n] > y have already been skipped above).
    if (n != 1)
    {
      for (uint32_t prime : filter_primes)
      {
        if (prime < n)
        {
          int64_t prime_index = prime_indexes[prime];
          bool expected = mu[n] != 0 && lpf[n] > prime;
          bool actual = factorTable[i] >
                        FactorTableD<uint16_t>::encode(prime_index);
          check(actual == expected);
        }
      }
    }

    not_coprime:;
  }

  // Compare single-threaded and parallel construction. The larger z
  // crosses FactorTableD's threading threshold.
  {
    int64_t parallel_y = 1000000;
    int64_t parallel_z = 5000001;
    FactorTableD<uint16_t> singleThread(parallel_y, parallel_z, 1);
    FactorTableD<uint16_t> multiThread(parallel_y, parallel_z, 2);
    int64_t size = singleThread.to_index(parallel_z) + 1;
    bool equal = true;

    for (int64_t i = 0; equal && i < size; i++)
    {
      equal = singleThread[i] == multiThread[i];

      if (equal && singleThread[i] != 0)
        equal = singleThread.mu(i) == multiThread.mu(i);
    }

    check(equal);
  }

  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
