///
/// @file   factor_table.cpp
/// @brief  FactorTable is a compressed lookup table
///         of mu (moebius) and lpf (least prime factor).
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <FactorTable.hpp>
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

  auto max = dist(gen);
  auto threads = max % 4;
  auto lpf = generate_lpf(max);
  auto mu = generate_moebius(max);

  FactorTable<int> factor_table(max, threads);

  for (int i = 2; i <= max; i++)
  {
    if (i % 2 != 0 &&
        i % 3 != 0 &&
        i % 5 != 0 &&
        i % 7 != 0)
    {
      if (mu[i] != 0)
      {
        cout << "mu(" << i << ") = " << mu[i];
        check(mu[i] == factor_table.mu(factor_table.get_index(i)));

        cout << "lpf(" << i << ") = " << lpf[i];
        check(lpf[i] <= factor_table.lpf(factor_table.get_index(i)) + (mu[i] == 1));
      }
    }
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
