///
/// @file  pi_meissel.cpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>

using namespace std;

namespace primecount {

/// Count the number of primes <= x using Meissel's formula.
/// Run time: O(x/(log x)^3)
/// Memory usage: O(x^0.5 / log x)
///
int64_t pi_meissel(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x);
  int64_t a = pi_legendre(y, 1);

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_meissel(x) ===" << endl;
    cout << "pi(x) = phi(x, a) + a - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "a = " << a << endl;
    cout << "threads = " << threads << endl;
  }

  int64_t p1 = phi(x, a, threads);
  int64_t p2 = P2(x, y, threads);
  int64_t sum = p1 + a - 1 - p2;

  return sum;
}

} // namespace
