///
/// @file   phi_tiny.cpp
/// @brief  Test the partial sieve function phi_tiny(x, a)
///         which counts the numbers <= x that are not divisible
///         by any of the first a primes with a <= 7.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <generate.hpp>

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

  int64_t max_a = PhiTiny::max_a();
  int64_t size = dist(gen);
  int64_t x = size - 1;

  auto primes = generate_n_primes<int>(max_a);
  vector<char> sieve(size, 1);

  for (int a = 1; a <= max_a; a++)
  {
    // remove primes[a] and its multiples
    for (int j = primes[a]; j <= x; j += primes[a])
      sieve[j] = 0;

    cout << "phi_tiny(" << x << ", " << a << ") = " << phi_tiny(x, a);
    check(phi_tiny(x, a) == count(sieve));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
