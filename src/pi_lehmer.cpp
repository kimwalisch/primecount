///
/// @file  pi_lehmer.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <pmath.hpp>

#include <stdint.h>

namespace primecount {

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2)) space.
///
int64_t pi_lehmer(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t x14 = iroot<4>(x);
  int64_t a = pi_meissel(x14, /* threads = */ 1);
  int64_t sum = 0;

  sum += phi(x, a, threads) + a - 1;
  sum -= P2 (x, a, threads);
  sum -= P3 (x, a, threads);

  return sum;
}

/// This version is for testing only, it uses a different P2(x, a)
/// implementation than pi_lehmer(x).
///
int64_t pi_lehmer2(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t x14 = iroot<4>(x);
  int64_t a = pi_meissel(x14, /* threads = */ 1);
  int64_t sum = 0;

  sum += phi(x, a, threads) + a - 1;
  sum -= P2 (x, x14);
  sum -= P3 (x, a, threads);

  return sum;
}

} // namespace primecount
