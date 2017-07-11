///
/// @file   alpha.cpp
/// @brief  Test the alpha tuning factor.
///         y = alpha * x^(1/3)
///         By computing pi(x) using different alpha tuning
///         factors we can make sure that all array sizes
///         (and other bounds) are accurate.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

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

  for (int i = 0; i < 20; i++)
  {
    int64_t x = dist(gen);
    int64_t res1 = pi_meissel(x);

    for (double alpha = 1; alpha <= iroot<6>(x); alpha++)
    {
      set_alpha(alpha);
      int64_t res2 = pi_deleglise_rivat(x);

      cout << "pi_deleglise_rivat(" << x << ") = " << res2;
      check(res1 == res2);
    }
  }

  min = ipow(10, 7);
  max = min * 2;
  uniform_int_distribution<int64_t> dist_small(min, max);

  for (int i = 0; i < 20; i++)
  {
    int64_t x = dist_small(gen);
    int64_t res1 = pi_meissel(x);

    for (double alpha = 1; alpha <= iroot<6>(x); alpha++)
    {
      set_alpha(alpha);
      int64_t res2 = pi_deleglise_rivat(x);

      cout << "pi_deleglise_rivat(" << x << ") = " << res2;
      check(res1 == res2);
    }
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
