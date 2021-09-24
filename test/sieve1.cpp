///
/// @file   sieve1.cpp
/// @brief  Test primecount's highly optimized modulo 30 sieve
///         of Eratosthenes implementation, specifically
///         Sieve::cross_off() and Sieve::count(low, high).
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <generate.hpp>
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

  auto primes = generate_primes<int>(isqrt(high));
  auto segment_size = Sieve::get_segment_size(high - low);

  Sieve sieve(low, segment_size, primes.size());
  std::vector<int> sieve2(high, 1);
  sieve2[0] = 0;

  for (size_t i = 1; i < primes.size(); i++)
  {
    if (primes[i] <= 5)
      sieve.pre_sieve(primes, i, low, high);
    else
      sieve.cross_off(primes[i], i);

    for (int j = primes[i]; j < high; j += primes[i])
      sieve2[j] = 0;

    if (primes[i] >= 5)
    {
      int start = dist(gen) % high;
      int stop  = dist(gen) % high;

      if (start > stop)
        std::swap(start, stop);

      uint64_t count = 0;

      for (int j = start; j <= stop; j++)
        count += sieve2[j];

      std::cout << "sieve.count(" << start << ", " << stop << ") = " << sieve.count(start, stop);
      check(count == sieve.count(start, stop));
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
