///
/// @file   S2_hard_xy.cpp
/// @brief  Test the computation of the hard special leaves
///         S2_hard(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
///         Note: when we set y = x^(1/3) then there are no
///         trivial and no easy special leaves which allows
///         us to test only the hard special leaves.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
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
  uniform_int_distribution<int> dist(1, 10000000);

  // test small x
  for (int i = 1; i < 30000; i++)
  {
    int64_t x = i;
    int64_t y = iroot<3>(x);
    int64_t pi_y = pi_legendre(y);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_n_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    cout << "S2_hard(" << x << ", " << y << ") = " << s2;
    check(s2 == S2_hard(x, y, z, c, Ri(x), 1));
  }

  // test random x
  for (int i = 0; i < 500; i++)
  {
    int64_t x = dist(gen);
    int64_t y = iroot<3>(x);
    int64_t pi_y = pi_legendre(y);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_n_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    cout << "S2_hard(" << x << ", " << y << ") = " << s2;
    check(s2 == S2_hard(x, y, z, c, Ri(x), 1));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
