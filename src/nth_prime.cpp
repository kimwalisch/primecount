///
/// @file  nth_prime.cpp
/// @brief Find the nth prime
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>

#include <stdint.h>
#include <string>
#include <array>

using namespace std;

namespace {

// Number of primes < 2^63
const int64_t max_n = 216289611853439384ll;

// primes[1] = 2, primes[2] = 3, ...
const array<int, 10> primes = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

}

namespace primecount {

/// Find the nth prime using a combination of the
/// Deleglise-Rivat prime counting algorithm and the
/// segmented sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/2))
///
int64_t nth_prime(int64_t n, int threads)
{
  if (n < 2)
    return primes[1];
  if (n < 10)
    return primes[n];

  if (n > max_n)
  {
    string msg("nth_prime(n): n must be <= " + to_string(max_n));
    throw primecount_error(msg);
  }

  int64_t prime = 0;

  if (n < 100000)
    prime = primesieve::nth_prime(n, 0);
  else
  {
    int64_t prime_approx = Ri_inverse(n);
    int64_t count_approx = pi(prime_approx, threads);

    if (count_approx < n)
      prime = primesieve::nth_prime(n - count_approx, prime_approx);
    else // count_approx >= n
      prime = primesieve::nth_prime(n - count_approx - 1, prime_approx + 1);
  }

  return prime;
}

} // namespace
