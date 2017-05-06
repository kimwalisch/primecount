///
/// @file  pi_legendre.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>

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

  // disable printing for pi_legendre()
  bool is_print = print_status();
  set_print_status(false);

  int64_t a = pi_primesieve(isqrt(x), /* threads = */ 1);
  int64_t sum = phi(x, a, threads) + a - 1;

  set_print_status(is_print);
  return sum;
}

} // namespace
