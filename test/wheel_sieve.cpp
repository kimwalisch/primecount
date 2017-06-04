///
/// @file   wheel_sieve.cpp
/// @brief  Test the Wheel class which is used to skip
///         multiples of 2, 3, 5 and 7 in the sieve of
///         Eratosthenes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Wheel.hpp>
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

  int low = dist(gen);
  auto primes = generate_primes<int>(isqrt(low));
  Wheel wheel(primes, primes.size(), low);

  for (size_t i = 5; i < primes.size(); i++)
  {
    int j = low;

    while (j % primes[i] != 0 ||
           j % 2 == 0 ||
           j % 3 == 0 ||
           j % 5 == 0 ||
           j % 7 == 0)
    {
      j++;
    }

    auto multiple = wheel[i].next_multiple;
    cout << "wheel.multiple(" << low << ", " << primes[i] << ") = " << multiple;
    check(multiple == j);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
