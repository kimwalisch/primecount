///
/// @file  nth_prime.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>

#include <stdint.h>

namespace primecount {

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/2)) space.
///
int64_t nth_prime(int64_t n, int threads)
{
  if (n < 1)
    n = 1;

  int64_t prime_approx = 0;
  int64_t count_approx = 0;

  if (n > 100000)
  {
    // Li_inverse(n) < nth_prime(n) for 7 <= n <= ~ 10^316
    prime_approx = Li_inverse(n);
    count_approx = pi(prime_approx, threads);
  }

  primesieve::set_num_threads(threads);
  int64_t prime = primesieve::parallel_nth_prime(n - count_approx, prime_approx);

  return prime;
}

} // namespace primecount
