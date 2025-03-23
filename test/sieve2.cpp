///
/// @file   sieve2.cpp
/// @brief  Test the return value of Sieve::cross_off(prime)
///         which returns the number of multiples of prime
///         that have been crossed off for the first time in
///         the sieve array.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::size_t;
using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(1000000, 2000000);

  int low = 0;
  int high = dist(gen);
  int sqrt_high = isqrt(high);
  auto primes = generate_primes<int32_t>(sqrt_high);

  uint64_t segment_size = high - low;
  segment_size = Sieve::align_segment_size(segment_size);
  Sieve sieve(low, segment_size, primes.size());
  std::vector<int> sieve2(high, 1);
  sieve2[0] = 0;

  for (size_t i = 1; i < primes.size(); i++)
  {
    uint64_t cnt1 = 0;
    uint64_t cnt2 = 0;
    uint64_t total1 = 0;
    uint64_t total2 = 0;

    if (primes[i] <= 5)
      sieve.pre_sieve(primes, i, low, high);
    else
    {
      uint64_t prev_count = sieve.get_total_count();
      sieve.cross_off_count(primes[i], i);
      cnt1 = prev_count - sieve.get_total_count();
      total1 = sieve.count(high - 1);
    }

    for (int j = primes[i]; j < high; j += primes[i])
    {
      cnt2 += sieve2[j];
      sieve2[j] = 0;
    }

    for (int j = 0; j < high; j++)
      total2 += sieve2[j];

    if (primes[i] > 5)
    {
      std::cout << "sieve.cross_off_count(" << i << ", " << primes[i] << ") = " << cnt1;
      check(cnt1 == cnt2);

      std::cout << "sieve.count(" << high - 1 << ") = " << total1;
      check((total1 == total2) && (total2 == sieve.get_total_count()));
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
