///
/// @file  pi_meissel.cpp
/// @brief Count the number of primes <= x using Meissel's formula.
///        Meissel's formula is an improved version of Legendre's
///        formula, it adds the P2(x, a) term which is the 2nd partial
///        sieve function and uses y=x^(1/3) instead of x^(1/2).
///
///        Meissel's formula:
///        pi(x) = pi(y) + phi(x, a) - 1 - P2(x, a)
///        with y = x^1/3, a = pi(y)
///
///        Please note that I think that Meissel's algorithm uses
///        O(x^(1/3)) memory instead of O(x^(1/2) / log(x)) found in
///        many papers. The memory usage is dominated by the segment
///        size (of the segmented sieve of Eratosthenes) in its P2
///        formula, which is O(sqrt(x^(2/3))) = O(x^(1/3)). One does
///        not need to store the primes up to x^(1/2) in P2 (as can be
///        seen in primecount's P2.cpp).
///
///        However, our implementation uses O(x^(1/2)) memory instead
///        of O(x^(1/3)) because our phi(x, a) implementation uses a
///        large pi(x) lookup table of size x^(1/2) in order to
///        improve performance.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

/// Count the number of primes <= x using Meissel's formula.
/// Run time: O(x/(log x)^3)
/// Memory usage: O(x^(1/2))
///
int64_t pi_meissel(int64_t x,
                   int threads,
                   bool is_print)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x);
  int64_t a = pi_noprint(y, threads);

  if (is_print)
  {
    print("");
    print("=== pi_meissel(x) ===");
    print("pi(x) = phi(x, a) + a - 1 - P2");
    print("x", x);
    print("y", y);
    print("a", a);
    print("threads", threads);
  }

  int64_t phi_xa = phi(x, a, threads, is_print);
  int64_t p2 = P2(x, y, a, threads, is_print);
  int64_t sum = phi_xa + a - 1 - p2;

  return sum;
}

} // namespace
