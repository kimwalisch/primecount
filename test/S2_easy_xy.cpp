///
/// @file   S2_easy_xy.cpp
/// @brief  Test the computation of the easy special leaves
///         S2_easy(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <generate.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <exception>
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
  try
  {
  
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(1, 100000000);

  // test small x
  for (int i = 1; i < 100000; i++)
  {
    int64_t x = i;
    double alpha = get_alpha_deleglise_rivat(x);
    int64_t x13 = iroot<3>(x);
    int64_t y = (int64_t) (alpha * x13);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2_easy = 0;

    auto primes = generate_primes<int64_t>(y);

    PiTable pi(y);
    int64_t pi_sqrty = pi[isqrt(y)];
    int64_t pi_x13 = pi[x13];

    for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
    {
      int64_t min_trivial = min(x / (primes[b] * primes[b]), y);
      int64_t min_sparse = max(z / primes[b], primes[b]);
      int64_t l = pi[min_trivial];

      for (; primes[l] > min_sparse; l--)
        s2_easy += pi[x / (primes[b] * primes[l])] - b + 2;
    }

    cout << "S2_easy(" << x << ", " << y << ") = " << s2_easy;
    check(s2_easy == S2_easy(x, y, z, c, 1));
  }

  // test random x
  for (int i = 0; i < 10000; i++)
  {
    int64_t x = dist(gen);
    double alpha = get_alpha_deleglise_rivat(x);
    int64_t x13 = iroot<3>(x);
    int64_t y = (int64_t) (alpha * x13);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2_easy = 0;

    auto primes = generate_primes<int64_t>(y);

    PiTable pi(y);
    int64_t pi_sqrty = pi[isqrt(y)];
    int64_t pi_x13 = pi[x13];

    for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
    {
      int64_t min_trivial = min(x / (primes[b] * primes[b]), y);
      int64_t min_sparse = max(z / primes[b], primes[b]);
      int64_t l = pi[min_trivial];

      for (; primes[l] > min_sparse; l--)
        s2_easy += pi[x / (primes[b] * primes[l])] - b + 2;
    }

    cout << "S2_easy(" << x << ", " << y << ") = " << s2_easy;
    check(s2_easy == S2_easy(x, y, z, c, 1));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;
  
  }
  catch (std::exception& e)
  {
    std::cout << "exception caught: " << e.what() << std::endl;
    std::cerr << "exception caught: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
