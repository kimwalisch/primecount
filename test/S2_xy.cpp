///
/// @file   S2_xy.cpp
/// @brief  Test the computation of the special leaves
///         S2(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <S.hpp>

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
  uniform_int_distribution<int> dist(1, 10000000);

  // test small x
  for (int i = 1; i < 30000; i++)
  {
    int64_t x = i;
    double alpha = get_alpha_deleglise_rivat(x);
    int64_t y = (int64_t) (alpha * iroot<3>(x));
    int64_t pi_y = pi_simple(y, 1);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    cout << "S2(" << x << ", " << y << ") = " << s2;

    check(s2 == S2_trivial(x, y, z, c) +
                S2_easy(x, y, z, c, 1) +
                S2_hard(x, y, z, c, Ri(x), 1));
  }

  // test random x
  for (int i = 0; i < 500; i++)
  {
    int64_t x = dist(gen);
    double alpha = get_alpha_deleglise_rivat(x);
    int64_t y = (int64_t) (alpha * iroot<3>(x));
    int64_t pi_y = pi_simple(y, 1);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    cout << "S2(" << x << ", " << y << ") = " << s2;

    check(s2 == S2_trivial(x, y, z, c) +
                S2_easy(x, y, z, c, 1) +
                S2_hard(x, y, z, c, Ri(x), 1));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
