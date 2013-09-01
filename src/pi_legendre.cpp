///
/// @file  pi_legendre.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "imath.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^0.5) space.
///
int64_t pi_legendre(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t a = pi_primesieve(isqrt(x), /* threads = */ 1);
  int64_t sum = phi(x, a, threads) + a - 1;

  return sum;
}

} // namespace primecount
