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

using namespace std;
using namespace primecount;

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(1000000, 2000000);

  int low = 0;
  int high = dist(gen);
  auto primes = generate_primes<int>(isqrt(high));

  Sieve sieve(low, high, primes.size());
  vector<int> sieve2(high, 1);
  sieve2[0] = 0;

  for (size_t i = 1; i < primes.size(); i++)
  {
    if (primes[i] <= 5)
      sieve.pre_sieve(i, low, high);
    else
      sieve.cross_off(i, primes[i]);

    for (int j = primes[i]; j < high; j += primes[i])
      sieve2[j] = 0;

    if (primes[i] >= 5)
    {
      int start = dist(gen) % high;
      int stop  = dist(gen) % high;

      if (start > stop)
        swap(start, stop);

      uint64_t count = 0;

      for (int j = start; j <= stop; j++)
        count += sieve2[j];

      cout << "sieve.count(" << start << ", " << stop << ") = " << sieve.count(start, stop);
      check(count == sieve.count(start, stop));
    }
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
