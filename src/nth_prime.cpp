///
/// @file  nth_prime.cpp
/// @brief Find the nth prime
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primecount.hpp>
#include <primesieve.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <string>

namespace {

const std::string max_n = "216289611853439384";

// primes[1] = 2, primes[2] = 3, ...
const int primes[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

}

namespace primecount {

/// Find the nth prime using a combination of the Deleglise-Rivat
/// prime counting algorithm and the segmented sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2), space: O(x^(1/2)).
///
int64_t nth_prime(int64_t n, int threads)
{
  n = std::max((int64_t) 1, n);

  if (n < 10)
    return primes[n];

  if (n > to_maxint(max_n))
    throw primecount_error("nth_prime(n): n must be <= " + max_n);

  int64_t prime = 0;
  int64_t prime_approx = 0;
  int64_t count_approx = 0;
  primesieve::set_num_threads(threads);

  if (n < 100000)
    prime = primesieve::parallel_nth_prime(n, 0);
  else
  {
    // Formula due to Dana Jacobsen:
    // Nth prime ≈ Li^-1(n) + Li^-1(sqrt(n)) / 4
    prime_approx = Li_inverse(n) + Li_inverse(isqrt(n)) / 4;
    count_approx = pi(prime_approx, threads);

    if (count_approx < n)
      prime = primesieve::parallel_nth_prime(n - count_approx, prime_approx);
    else /* count_approx >= n */
      prime = primesieve::parallel_nth_prime(n - count_approx - 1, prime_approx + 1);
  }

  return prime;
}

} // namespace primecount
