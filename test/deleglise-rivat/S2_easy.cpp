///
/// @file   S2_easy.cpp
/// @brief  Test the computation of the easy special leaves
///         S2_easy(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
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

  threads = get_num_threads();

  {
    // Test S2_easy(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 178815;
    int64_t z = 55923720;
    int64_t c = 8;
    int64_t res1 = S2_easy(x, y, z, c, threads);
    int64_t res2 = 60888055472ll;

    std::cout << "S2_easy(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test S2_easy(1e14) and compare with known correct value
    int128_t x = 100000000000000ll;
    int64_t y = 494134;
    int64_t z = 202374254;
    int64_t c = 8;
    int128_t res1 = S2_easy(x, y, z, c, threads);
    int128_t res2 = 617442826127ll;

    std::cout << "S2_easy(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
