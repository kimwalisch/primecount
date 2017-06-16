///
/// @file   bit_sieve.cpp
/// @brief  Test the BitSieve class
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve.hpp>
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
  uniform_int_distribution<int> dist(500000, 1000000);

  int size = dist(gen);
  int low = 1;

  auto primes = generate_primes<int>(isqrt(size));
  BitSieve bitSieve(size);
  vector<int> sieve(size, 1);
  size_t pre_sieve = 7;

  for (int j = primes[1] - low; j < size; j += primes[1])
    sieve[j] = 0;

  for (size_t i = 2; i < primes.size(); i++)
  {
    if (i <= pre_sieve)
      bitSieve.pre_sieve(i, low);

    for (int j = primes[i] - low; j < size; j += primes[i] * 2)
    {
      sieve[j] = 0;
      if (i > pre_sieve)
        bitSieve.unset(j);
    }

    int start = dist(gen) % size;
    int stop  = dist(gen) % size;

    if (start > stop)
      swap(start, stop);

    uint64_t count = 0;

    for (int j = start; j <= stop; j++)
      count += sieve[j];

    cout << "bitSieve.count(" << start << ", " << stop << ") = " << bitSieve.count(start, stop);
    check(count == bitSieve.count(start, stop));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
