///
/// @file  pi_lehmer.cpp
/// @brief Count the number of primes <= x using Lehmer's formula.
///        Lehmer's formula is an improved version of Meissel's
///        formula, it adds the P3(x, a) term which is the 3rd partial
///        sieve function and uses y=x^(1/4) instead of x^(1/3).
///
///        Lehmer's formula:
///        pi(x) = pi(y) + phi(x, a) - 1 - P2(x, a) - P3(x, a)
///        with y = x^1/4, a = pi(y)
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

/// Count the number of primes <= x using Lehmer's formula.
/// Run time: O(x/(log x)^4)
/// Memory usage: O(x^(1/2))
///
int64_t pi_lehmer(int64_t x,
                  int threads,
                  bool is_print)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<4>(x);
  int64_t a = pi_noprint(y, threads);

  if (is_print)
  {
    print("");
    print("=== pi_lehmer(x) ===");
    print("pi(x) = phi(x, a) + a - 1 - P2 - P3");
    print("x", x);
    print("y", y);
    print("a", a);
    print("threads", threads);
  }

  int64_t phi_xa = phi(x, a, threads, is_print);
  int64_t p2 = P2(x, y, threads, is_print);
  int64_t p3 = P3(x, y, threads, is_print);
  int64_t sum = phi_xa + a - 1 - p2 - p3;

  return sum;
}

} // namespace
