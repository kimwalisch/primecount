///
/// @file  nth_prime.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "isqrt.h"

#include <primecount.h>
#include <primesieve/soe/ParallelPrimeSieve.h>
#include <stdint.h>

namespace primecount {

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x/(log x)^4) operations, O(x^0.5) space.
///
int64_t nth_prime(int64_t n, int threads)
{
  if (n < 1)
    return 0;

  int64_t prime_approx = 0;
  int64_t count_approx = 0;

  if (n > 100000)
  {
    // Li_inverse(n) < nth_prime(n) for 7 <= n <= ~ 10^316
    prime_approx = Li_inverse(n);
    count_approx = pi(prime_approx, threads);
  }

  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS)
    pps.setNumThreads( threads );
  int64_t prime = pps.nthPrime(prime_approx + 1, n - count_approx);

  return prime;
}

} // namespace primecount
