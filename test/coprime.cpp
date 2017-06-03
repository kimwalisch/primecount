///
/// @file   coprime.cpp
/// @brief  Test the partial sieve function phi(x, a)
///         which counts the numbers <= x that are not divisible
///         by any of the first a primes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
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

// Count the number of unsieved elements
int count(vector<char>& sieve)
{
  int cnt = 0;

  for (size_t i = 1; i < sieve.size(); i++)
    cnt += sieve[i];

  return cnt;
}

int main()
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(10000000, 20000000);

  int size = dist(gen);
  int x = size - 1;

  auto primes = generate_primes<int>(isqrt(x));
  vector<char> sieve(size, 1);

  for (size_t a = 1; a < primes.size(); a++)
  {
    // remove primes[a] and its multiples
    for (int j = primes[a]; j <= x; j += primes[a])
      sieve[j] = 0;

    cout << "phi(" << x << ", " << a << ") = " << phi(x, a);
    check(phi(x, a) == count(sieve));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
