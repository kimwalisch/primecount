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
#include <PiTable.hpp>
#include <primecount-config.hpp>

#include <stdint.h>
#include <algorithm>
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
  check(FactorTable<uint16_t>::to_index(FactorTable<uint16_t>::max()) > uint32_max);
  check(FactorTable<uint16_t>::encode(32766) == uint16_max - 2);
  check(FactorTable<uint16_t>::encode(32767) == uint16_max);

  bool rejected = false;
  try
  {
    PiTable pi(100, 1);
    auto primes = pi.get_primes<uint32_t>(100, 1);
    FactorTable<uint16_t>(FactorTable<uint16_t>::max() + 1, pi, primes, 1);
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
  auto threads = std::min(max % 4, get_num_threads());
  PiTable pi(max, threads);
  auto primes = pi.get_primes<uint32_t>(max, threads);
  auto lpf = generate_lpf(max, primes);
  auto mu = generate_moebius(max, primes);
  Vector<uint32_t> prime_indexes(max + 1);

  for (std::size_t i = 1; i < primes.size(); i++)
    prime_indexes[primes[i]] = (uint32_t) i;

  // FactorTable<uint16_t>::max() requires p(32767) = 386083
  check(primes[32767] == 386083);

  auto primes64 = pi.get_primes<int64_t>(max, threads);
  FactorTable<uint16_t> factorTable16(max, pi, primes, threads);
  FactorTable<uint32_t> factorTable32(max, pi, primes64, threads);
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

  // Compare single-threaded and parallel construction
  {
    // Numbers per L2_CACHE_SIZE sized segment of a
    // FactorTable<uint16_t>, FactorTable<uint32_t> segments
    // contain half as many numbers. Hence a single thread
    // initializes both factor tables in multiple segments.
    int64_t segment_size = L2_CACHE_SIZE / sizeof(uint16_t) * 2310ull / 480;
    int64_t parallel_y = segment_size * 2 + 1;
    int parallel_threads = std::min(2, get_num_threads());

    // FactorTable requires the primes <= y / 13
    PiTable parallel_pi(parallel_y / 13, 1);
    auto parallel_primes = parallel_pi.get_primes<uint32_t>(parallel_y / 13, 1);
    FactorTable<uint16_t> singleThread16(parallel_y, parallel_pi, parallel_primes, 1);
    FactorTable<uint16_t> multiThread16(parallel_y, parallel_pi, parallel_primes, parallel_threads);
    FactorTable<uint32_t> singleThread32(parallel_y, parallel_pi, parallel_primes, 1);
    FactorTable<uint32_t> multiThread32(parallel_y, parallel_pi, parallel_primes, parallel_threads);
    int64_t size = singleThread16.to_index(parallel_y) + 1;
    bool equal = true;

    for (int64_t i = 0; equal && i < size; i++)
    {
      equal = singleThread16[i] == multiThread16[i];
      equal = equal && singleThread32[i] == multiThread32[i];

      if (equal && singleThread16[i] != 0)
        equal = singleThread16.mu(i) == multiThread16.mu(i);
    }

    check(equal);

    // Check the single-threaded factor tables against the
    // brute-force lpf(n) and mu(n) values. The numbers <= max_n
    // span the first segment boundary of both factor tables.
    int64_t max_n = segment_size + segment_size / 32;
    PiTable pi_y(max_n, 1);
    auto primes_y = pi_y.get_primes<uint32_t>(max_n, 1);
    auto lpf_y = generate_lpf(max_n, primes_y);
    auto mu_y = generate_moebius(max_n, primes_y);
    int64_t max_idx = singleThread16.to_index(max_n);

    for (int64_t i = 0; i <= max_idx; i++)
    {
      int64_t n = singleThread16.to_number(i);
      bool OK;

      if (n == 1)
        OK = (singleThread16[i] == uint16_max - 1 &&
              singleThread32[i] == uint32_max - 1);
      else if (lpf_y[n] == n)
        OK = (singleThread16[i] == uint16_max &&
              singleThread32[i] == uint32_max);
      else if (mu_y[n] == 0)
        OK = (singleThread16[i] == 0 && singleThread32[i] == 0);
      else
      {
        int64_t expected = pi_y[lpf_y[n]] * 2 + (mu_y[n] == -1);
        OK = (singleThread16[i] == expected && singleThread32[i] == expected);
      }

      if (!OK)
        std::cout << "lpf(" << n << ") = " << lpf_y[n];

      check(OK);
    }
  }

  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
