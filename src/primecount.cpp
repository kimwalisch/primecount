///
/// @file   primecount.cpp
/// @brief  primecount C++ API
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <utils.hpp>

#include <stdint.h>

namespace {

// By default primecount uses all CPU cores.
int threads = primecount::MAX_THREADS;

}

namespace primecount {

// Set the number of threads.
void set_num_threads(int num_threads)
{
  threads = validate_threads(num_threads);
}

// Get the currently set number of threads.
int get_num_threads()
{
  return validate_threads(threads);
}

/// This is an alias for the fastest prime counting implementation
/// within primecount.
///
int64_t pi(int64_t x)
{
  return pi(x, threads);
}

/// This is an alias for the fastest prime counting implementation
/// within primecount.
///
int64_t pi(int64_t x, int threads)
{
  return pi_deleglise_rivat(x, threads);
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log x) space.
///
int64_t pi_deleglise_rivat(int64_t x)
{
  return pi_deleglise_rivat(x, threads);
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log x) space.
///
int64_t pi_deleglise_rivat(int64_t x, int threads)
{
  return pi_deleglise_rivat_parallel3(x, threads);
}

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^(1/2)) space.
///
int64_t pi_legendre(int64_t x)
{
  return pi_legendre(x, threads);
}

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2)) space.
///
int64_t pi_lehmer(int64_t x)
{
  return pi_lehmer(x, threads);
}

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo(int64_t x)
{
  return pi_lmo(x, threads);
}

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo(int64_t x, int threads)
{
  return pi_lmo_parallel3(x, threads);
}

/// Calculate the number of primes below x using Meissel's formula.
/// Run time: O(x/(log x)^3) operations, O(x^(1/2) / log x) space.
///
int64_t pi_meissel(int64_t x)
{
  return pi_meissel(x, threads);
}

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x)
{
  return pi_primesieve(x, threads);
}

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t nth_prime(int64_t n)
{
  return nth_prime(n, threads);
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a)
{
  return phi(x, a, threads);
}

} // namespace primecount
