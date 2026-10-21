///
/// @file   FactorTable.cpp
/// @brief  FactorTable is a compressed lookup table of mu
///         (moebius) and lpf (least prime factor).
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

// factorTable.mu(n) = 0 is disabled by default for performance
// reasons, we only enable it for testing.
#define ENABLE_MU_0_TESTING

#include <FactorTable.hpp>
#include <generate_primes.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
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
  int64_t uint16_max = pstd::numeric_limits<uint16_t>::max();
  int64_t uint32_max = pstd::numeric_limits<uint32_t>::max();
  check(FactorTable<uint16_t>::max() == 149060082888);
  check(FactorTable<uint32_t>::max() == pstd::numeric_limits<int64_t>::max());
  check(FactorTable<uint16_t>::to_index(
          FactorTable<uint16_t>::max()) > UINT32_MAX);
  check(FactorTable<uint16_t>::encode(32766) == uint16_max - 2);
  check(FactorTable<uint16_t>::encode(32767) == uint16_max);

  bool rejected = false;
  try
  {
    FactorTable<uint16_t>(FactorTable<uint16_t>::max() + 1, 1);
  }
  catch (const primecount_error&)
  {
    rejected = true;
  }
  check(rejected);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(500000, 1000000);

  auto max = dist(gen);
  auto threads = max % 4;
  auto primes = generate_primes<uint32_t>(max);
  auto lpf = generate_lpf(max, primes);
  auto mu = generate_moebius(max, primes);
  Vector<uint32_t> prime_indexes(max + 1);

  for (std::size_t i = 1; i < primes.size(); i++)
    prime_indexes[primes[i]] = (uint32_t) i;

  // FactorTable<uint16_t>::max() requires p(32767) = 386083
  check(primes[32767] == 386083);

  FactorTable<uint16_t> factorTable16(max, threads);
  FactorTable<uint32_t> factorTable32(max, threads);
  int64_t limit = factorTable16.first_coprime();
  Array<uint32_t, 4> filter_primes = { 13, 17, 101, 9973 };

  for (int64_t n = 1; n <= max; n++)
  {
    int64_t i = factorTable16.to_index(n);
    bool is_prime = (lpf[n] == n);
    bool mu_OK;
    bool lpf_OK;

    // Check if n is coprime to the primes < limit
    for (int64_t j = 1; primes[j] < limit; j++)
      if (n % primes[j] == 0)
        goto not_coprime;

    mu_OK = (mu[n] == factorTable16.mu(i) &&
             mu[n] == factorTable32.mu(i));

    if (!mu_OK)
      std::cout << "mu(" << n << ") = " << factorTable16.mu(i)
                << ", " << factorTable32.mu(i);

    check(mu_OK);

    // factorTable[index] is a combination of the mu(n) (Möbius
    // function) and lpf(n) (least prime factor) functions.
    // It returns (with n = to_number(index)):
    //
    // 1) INT_MAX - 1    if n = 1
    // 2) INT_MAX        if n is a prime
    // 3) 0              if moebius(n) = 0
    // 4) 2*pi(lpf)      if moebius(n) = 1
    // 5) 2*pi(lpf) + 1  if moebius(n) = -1

    if (n == 1)
      lpf_OK = (factorTable16[i] == uint16_max - 1 &&
                factorTable32[i] == uint32_max - 1);
    else if (is_prime)
      lpf_OK = (factorTable16[i] == uint16_max &&
                factorTable32[i] == uint32_max);
    else if (mu[n] == 0)
      lpf_OK = (factorTable16[i] == 0 && factorTable32[i] == 0);
    else
    {
      int64_t expected = prime_indexes[lpf[n]] * 2 + (mu[n] == -1);
      lpf_OK = (factorTable16[i] == expected && factorTable32[i] == expected);
    }

    if (!lpf_OK)
      std::cout << "lpf(" << n << ") = " << lpf[n];

    check(lpf_OK);

    // Check that "factor[n] > encode(PrimePi(prime))" is equivalent
    // to "mu[n] != 0 && lpf[n] > prime" for a few prime thresholds.
    if (n != 1)
    {
      for (uint32_t prime : filter_primes)
      {
        if (prime < n)
        {
          int64_t prime_index = prime_indexes[prime];
          bool expected = mu[n] != 0 && lpf[n] > prime;
          check((factorTable16[i] > FactorTable<uint16_t>::encode(prime_index)) == expected);
          check((factorTable32[i] > FactorTable<uint32_t>::encode(prime_index)) == expected);
        }
      }
    }

    not_coprime:;
  }

  // Compare single-threaded and parallel construction. The larger y
  // crosses FactorTable's threading threshold.
  {
    int64_t parallel_y = 5000001;
    FactorTable<uint16_t> singleThread(parallel_y, 1);
    FactorTable<uint16_t> multiThread(parallel_y, 2);
    FactorTable<uint32_t> singleThread32(parallel_y, 1);
    FactorTable<uint32_t> multiThread32(parallel_y, 2);
    int64_t size = singleThread.to_index(parallel_y) + 1;
    bool equal = true;

    for (int64_t i = 0; equal && i < size; i++)
    {
      equal = singleThread[i] == multiThread[i];
      equal = equal && singleThread32[i] == multiThread32[i];

      if (equal && singleThread[i] != 0)
        equal = singleThread.mu(i) == multiThread.mu(i);
    }

    check(equal);
  }

  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
