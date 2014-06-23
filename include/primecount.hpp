///
/// @file  primecount.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_HPP
#define PRIMECOUNT_HPP

#include <stdint.h>

#define PRIMECOUNT_VERSION "0.21"
#define PRIMECOUNT_VERSION_MAJOR 0
#define PRIMECOUNT_VERSION_MINOR 21

namespace primecount {

/// This is an alias for the fastest prime counting implementation
/// within primecount.
///
int64_t pi(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat(int64_t x);

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^(1/2)) space.
///
int64_t pi_legendre(int64_t x);

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2)) space.
///
int64_t pi_lehmer(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo(int64_t x);

/// Calculate the number of primes below x using Meissel's formula.
/// Run time: O(x/(log x)^3) operations, O(x^(1/2) / log x) space.
///
int64_t pi_meissel(int64_t x);

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x);

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t nth_prime(int64_t n);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a);

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes below x.
/// @post Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
int64_t Li(int64_t x);

/// Calculate the inverse logarithmic integral Li^-1(x) which is
/// a very accurate approximation of the nth prime.
/// @post Li_inverse(x) < nth_prime(x) for 7 <= x <= ~ 10^316
///
int64_t Li_inverse(int64_t x);

// Set the number of threads.
void set_num_threads(int num_threads);

// Get the currently set number of threads.
int get_num_threads();

/// Test all prime counting function implementations.
/// @return true if success else false.
///
bool test();

} // namespace primecount

#endif
