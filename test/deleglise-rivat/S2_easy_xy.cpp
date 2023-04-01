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

#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <S.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::min;
using std::max;
using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  int64_t max_x = 100000;
  double max_alpha = get_alpha_deleglise_rivat(max_x);
  int64_t max_x13 = iroot<3>(max_x);
  int64_t max_y = (int64_t) (max_alpha * max_x13);
  int threads = 1;

  {
    auto primes = generate_primes<int64_t>(max_y);
    PiTable pi(max_y, threads);

    // test small x
    for (int64_t i = 1; i < max_x; i++)
    {
      int64_t x = i;
      double alpha = get_alpha_deleglise_rivat(x);
      int64_t x13 = iroot<3>(x);
      int64_t y = (int64_t) (alpha * x13);
      int64_t z = x / y;
      int64_t c = PhiTiny::get_c(y);
      int64_t pi_sqrty = pi[isqrt(y)];
      int64_t pi_x13 = pi[x13];
      int64_t s2_easy = 0;

      for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
      {
        int64_t min_trivial = min(x / (primes[b] * primes[b]), y);
        int64_t min_sparse = max(z / primes[b], primes[b]);
        int64_t l = pi[min_trivial];

        for (; primes[l] > min_sparse; l--)
          s2_easy += pi[x / (primes[b] * primes[l])] - b + 2;
      }

      std::cout << "S2_easy(" << x << ", " << y << ") = " << s2_easy;
      check(s2_easy == S2_easy(x, y, z, c, threads));
    }
  }

  max_x = 100000000;
  max_alpha = get_alpha_deleglise_rivat(max_x);
  max_x13 = iroot<3>(max_x);
  max_y = (int64_t) (max_alpha * max_x13);

  {
    auto primes = generate_primes<int64_t>(max_y);
    PiTable pi(max_y, threads);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(1, max_x);

    // test random x
    for (int64_t i = 0; i < 10000; i++)
    {
      int64_t x = dist(gen);
      double alpha = get_alpha_deleglise_rivat(x);
      int64_t x13 = iroot<3>(x);
      int64_t y = (int64_t) (alpha * x13);
      int64_t z = x / y;
      int64_t c = PhiTiny::get_c(y);
      int64_t pi_sqrty = pi[isqrt(y)];
      int64_t pi_x13 = pi[x13];
      int64_t s2_easy = 0;

      for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
      {
        int64_t min_trivial = min(x / (primes[b] * primes[b]), y);
        int64_t min_sparse = max(z / primes[b], primes[b]);
        int64_t l = pi[min_trivial];

        for (; primes[l] > min_sparse; l--)
          s2_easy += pi[x / (primes[b] * primes[l])] - b + 2;
      }

      std::cout << "S2_easy(" << x << ", " << y << ") = " << s2_easy;
      check(s2_easy == S2_easy(x, y, z, c, threads));
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
