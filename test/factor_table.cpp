///
/// @file   factor_table.cpp
/// @brief  FactorTable is a compressed lookup table
///         of mu (moebius) and lpf (least prime factor).
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <FactorTable.hpp>
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

int main()
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(500000, 1000000);

  auto max = dist(gen);
  auto threads = max % 4;
  auto lpf = generate_lpf(max);
  auto mu = generate_moebius(max);

  vector<int> small_primes = { 2, 3, 5, 7, 11, 13, 17, 19 };
  FactorTable<int> factor_table(max, threads);
  int64_t limit = factor_table.get_first_coprime();

  for (int n = 1; n <= max; n++)
  {
    // Check if n is coprime to the primes < limit
    for (int p : small_primes)
    {
      if (p >= limit)
        break;
      if (n % p == 0)
        goto not_coprime;
    }

    // For performance reasons FactorTable does not support (mu[n] == 0)
    if (mu[n] != 0)
    {
      int64_t i = factor_table.get_index(n);

      cout << "mu(" << n << ") = " << mu[n];
      check(mu[n] == factor_table.mu(i));

      cout << "lpf(" << n << ") = " << lpf[n];
      check(lpf[n] <= factor_table.lpf(i) + (mu[n] == 1));
    }

    not_coprime:;
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
