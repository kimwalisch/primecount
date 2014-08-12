///
/// @file  nth_prime.cpp
/// @brief Find the nth prime
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <pmath.hpp>

#include <stdint.h>

namespace primecount {

/// Find the nth prime using a combination of the Deleglise-Rivat
/// prime counting algorithm and the segmented sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2), space: O(x^(1/2)).
///
int64_t nth_prime(int64_t n, int threads)
{
  if (n < 1)
    n = 1;

  int64_t prime = 0;
  int64_t prime_approx = 0;
  int64_t count_approx = 0;

  primesieve::set_num_threads(threads);

  if (n < 100000)
    prime = primesieve::parallel_nth_prime(n, 0);
  else
  {
    // Formula due to Dana Jacobsen:
    // Nth prime ~ Li^-1(n) + Li^-1(sqrt(n)) / 4
    prime_approx = Li_inverse(n) + Li_inverse(isqrt(n)) / 4;
    count_approx = pi(prime_approx, threads);

    if (count_approx <= n)
      prime = primesieve::parallel_nth_prime(n - count_approx, prime_approx);
    else /* count_approx > n */
      prime = primesieve::parallel_nth_prime(n - count_approx - 1, prime_approx + 1);
  }

  return prime;
}

} // namespace primecount
