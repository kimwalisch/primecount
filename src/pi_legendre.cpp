///
/// @file  pi_legendre.cpp
/// @brief Count the number of primes <= x using Legendre's formula.
///        Legendre's prime counting algorithm is the simplest
///        combinatorial algorithm for counting the number of primes.
///        All other formulas (e.g. Meissel's formula, Lehmer's
///        formula, ...) are extensions of Legendre's formula that
///        run faster but are also more complex.
///
///        Legendre's formula:
///        pi(x) = pi(y) + phi(x, pi(y)) - 1
///        with y = x^(1/2)
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <print.hpp>
#include <isqrt.hpp>

#include <stdint.h>

namespace primecount {

/// Count the number of primes <= x using Legendre's formula.
/// Run time: O(x)
/// Memory usage: O(x^(1/2))
///
int64_t pi_legendre(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t a = pi_legendre(isqrt(x), 1);
  int64_t sum = phi(x, a, threads) + a - 1;

  return sum;
}

/// The default pi_legendre(x) implementation does not support
/// printing as it is used all over the place in many other
/// formulas (e.g. for initializing the size of arrays) and we
/// don't want to print info about that.
///
int64_t pi_legendre_print(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t a = pi_legendre(isqrt(x), 1);

  print("");
  print("=== pi_legendre(x) ===");
  print("pi(x) = phi(x, a) + a - 1");
  print("x", x);
  print("a", a);
  print("threads", threads);

  int64_t sum = phi_print(x, a, threads) + a - 1;

  return sum;
}

} // namespace
