///
/// @file  primecount.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_H
#define PRIMECOUNT_H

#include <stdint.h>

#define PRIMECOUNT_VERSION "0.5"
#define PRIMECOUNT_MAJOR_VERSION 0
#define PRIMECOUNT_MINOR_VERSION 5
#define PRIMECOUNT_YEAR "2013"

namespace primecount {

enum {
  /// Uses all CPU cores, do not modify.
  MAX_THREADS = -1
};

/// This is an alias for the fastest prime counting implementation
/// within primecount.
///
int64_t pi(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^0.5) space.
///
int64_t pi_primesieve(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^0.5) space.
///
int64_t pi_legendre(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using Meissel's formula.
/// Run time: O(x/(log x)^3) operations, O(x^0.5) space.
///
int64_t pi_meissel(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^0.5) space.
///
int64_t pi_lehmer(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
///
int64_t pi_lmo_simple(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
///
int64_t pi_lmo(int64_t x, int threads = MAX_THREADS);

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x/(log x)^4) operations, O(x^0.5) space.
///
int64_t nth_prime(int64_t n, int threads = MAX_THREADS);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads = MAX_THREADS);

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes below x.
/// @post Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
int64_t Li(int64_t);

/// Calculate the inverse logarithmic integral Li^-1(x) which is
/// a very accurate approximation of the nth prime.
/// @post Li_inverse(x) < nth_prime(x) for 7 <= x <= ~ 10^316
///
int64_t Li_inverse(int64_t);

} // namespace primecount

#endif
