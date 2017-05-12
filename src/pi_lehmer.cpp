///
/// @file  pi_lehmer.cpp
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

using namespace std;

namespace primecount {

/// Count the number of primes <= x using Lehmer's formula.
/// Run time: O(x/(log x)^4)
/// Memory usage: O(x^(1/2))
///
int64_t pi_lehmer(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<4>(x);
  int64_t a = pi_legendre(y);

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_lehmer(x) ===" << endl;
    cout << "pi(x) = phi(x, a) + a - 1 - P2 - P3" << endl;
    cout << "x = " << x << endl;
    cout << "a = " << a << endl;
    cout << "threads = " << threads << endl;
  }

  int64_t p1 = phi(x, a, threads);
  int64_t p2 = P2(x, y, threads);
  int64_t p3 = P3(x, a, threads);
  int64_t sum = p1 + a - 1 - p2 - p3;

  return sum;
}

} // namespace
