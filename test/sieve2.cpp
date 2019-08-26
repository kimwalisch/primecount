///
/// @file   sieve2.cpp
/// @brief  Test the return value of Sieve::cross_off(prime)
///         which returns the number of multiples of prime
///         that have been crossed off for the first time in
///         the sieve array.
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
  auto segment_size = Sieve::get_segment_size(high - low);

  Sieve sieve(low, segment_size, primes.size());
  vector<int> sieve2(high, 1);
  sieve2[0] = 0;

  for (size_t i = 1; i < primes.size(); i++)
  {
    int64_t cnt1 = 0;
    int64_t cnt2 = 0;

    if (primes[i] <= 5)
      sieve.pre_sieve(primes, i, low, high);
    else
      cnt1 = sieve.cross_off_count(primes[i], i);

    for (int j = primes[i]; j < high; j += primes[i])
    {
      cnt2 += sieve2[j];
      sieve2[j] = 0;
    }

    if (primes[i] > 5)
    {
      cout << "sieve.cross_off_count(" << i << ", " << primes[i] << ") = " << cnt1;
      check(cnt1 == cnt2);
    }
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
