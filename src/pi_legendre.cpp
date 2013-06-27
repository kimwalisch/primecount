///
/// @file  pi_legendre.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "isqrt.h"
#include "phi.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_legendre(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t a = pi_primesieve(isqrt(x), threads);
  int64_t sum = a + phi(x, a, threads) - 1;

  return sum;
}

} // namespace primecount
