///
/// @file   alpha_z.cpp
/// @brief  Test the alpha_z tuning factor in Gourdon's algorithm.
///         z = y * alpha_z
///         By computing pi(x) using different alpha_z tuning
///         factors we can make sure that all array sizes
///         (and other bounds) are accurate.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
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

  int64_t min = ipow(10, 9);
  int64_t max = min * 2;
  uniform_int_distribution<int64_t> dist(min, max);
  int threads = get_num_threads();

  for (int i = 0; i < 20; i++)
  {
    int64_t x = dist(gen);
    int64_t res1 = pi_meissel(x, threads);

    for (double alpha_z = 1; alpha_z <= iroot<6>(x); alpha_z++)
    {
      set_alpha_z(alpha_z);
      int64_t res2 = pi_gourdon(x, threads);

      cout << "pi_gourdon(" << x << ") = " << res2;
      check(res1 == res2);
    }
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
