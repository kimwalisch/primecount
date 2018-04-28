///
/// @file  pi_meissel.cpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

/// Count the number of primes <= x using Meissel's formula.
/// Run time: O(x/(log x)^3)
/// Memory usage: O(x^(1/2))
///
int64_t pi_meissel(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x);
  int64_t a = pi_legendre(y);

  print("");
  print("=== pi_meissel(x) ===");
  print("pi(x) = phi(x, a) + a - 1 - P2");
  print("x", x);
  print("y", y);
  print("a", a);
  print("threads", threads);

  double time = get_wtime();
  int64_t p1 = phi(x, a, threads);

  print("");
  print("=== phi(x, a) ===");
  print("Count the numbers <= x coprime to the first a primes");
  print("phi", p1, time);

  int64_t p2 = P2(x, y, threads);
  int64_t sum = p1 + a - 1 - p2;

  return sum;
}

} // namespace
