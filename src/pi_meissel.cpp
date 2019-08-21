///
/// @file  pi_meissel.cpp
/// @brief Count the number of primes <= x using Meissel's formula.
///        Meissel's formula is an improvement over Legendre's formula,
///        it adds the P2(x, a) term which is the 2nd partial sieve
///        function and sets y=x^(1/3) instead of x^(1/2).
///
///        Meissel's formula:
///        pi(x) = pi(y) + phi(x, a) - 1 - P2(x, a)
///        with y = x^1/3, a = pi(y)
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
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
  int64_t a = pi_legendre(y, threads);

  print("");
  print("=== pi_meissel(x) ===");
  print("pi(x) = phi(x, a) + a - 1 - P2");
  print("x", x);
  print("y", y);
  print("a", a);
  print("threads", threads);

  print("");
  print("=== phi(x, a) ===");
  print("Count the numbers <= x coprime to the first a primes");
  double time = get_time();
  int64_t p1 = phi(x, a, threads);
  print("phi", p1, time);

  int64_t p2 = P2(x, y, threads);
  int64_t sum = p1 + a - 1 - p2;

  return sum;
}

} // namespace
