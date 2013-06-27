///
/// @file  pi_lehmer.cpp
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

int64_t pi_lehmer(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t a = pi_legendre(isqrt4(x));
  int64_t b = pi_legendre(isqrt(x));
  int64_t c = pi_legendre(isqrt3(x));

  int64_t sum = 0;

  sum += phi(x, a, threads);
  sum += phi2(x, a, b, isqrt(x), threads);
  sum += phi3(x, a, c, isqrt(x), threads);

  return sum;
}

} // namespace primecount
