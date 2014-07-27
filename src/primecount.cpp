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
#include <calculator.hpp>
#include <ptypes.hpp>
#include <utils.hpp>

#include <sstream>
#include <string>
#include <stdint.h>

namespace {

// By default primecount uses all CPU cores.
int threads_ = primecount::MAX_THREADS;

}

namespace primecount {

// Set the number of threads.
void set_num_threads(int threads)
{
  threads_ = validate_threads(threads);
}

// Get the currently set number of threads.
int get_num_threads()
{
  return validate_threads(threads_);
}

/// Alias for the fastest prime counting function in primecount.
int64_t pi(int64_t x)
{
  return pi(x, threads_);
}

/// Alias for the fastest prime counting function in primecount.
int64_t pi(int64_t x, int threads)
{
  return pi_deleglise_rivat(x, threads);
}

#ifdef HAVE_INT128_T

/// Alias for the fastest prime counting function in primecount.
int128_t pi(int128_t x)
{
  // TODO
  return -1;
}

/// Alias for the fastest prime counting function in primecount.
int128_t pi(int128_t x, int threads)
{
  // TODO
  return -1;
}

#endif /* HAVE_INT128_T */

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like 10^12.
/// @pre   x  <= primecount::max().
///
std::string pi(const std::string& x)
{
  return pi(x, threads_);
}

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like 10^12.
/// @pre   x  <= primecount::max().
///
std::string pi(const std::string& x, int threads)
{
  maxint_t n = calculator::eval<maxint_t>(x);
  maxint_t pin = pi(n, threads);
  std::ostringstream oss;
  oss << pin;
  return oss.str();
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat(int64_t x)
{
  return pi_deleglise_rivat(x, threads_);
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
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
  return pi_legendre(x, threads_);
}

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2)) space.
///
int64_t pi_lehmer(int64_t x)
{
  return pi_lehmer(x, threads_);
}

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo(int64_t x)
{
  return pi_lmo(x, threads_);
}

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
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
  return pi_meissel(x, threads_);
}

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x)
{
  return pi_primesieve(x, threads_);
}

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/2)) space.
///
int64_t nth_prime(int64_t n)
{
  return nth_prime(n, threads_);
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a)
{
  return phi(x, a, threads_);
}

/// Returns the largest integer that can be used with
/// pi(std::string x). The return type is a string as max may be a
/// 128-bit integer which is not supported by all compilers.
///
std::string max()
{
#ifdef HAVE_INT128_T
  // TODO
  return "-1";
#else
  std::ostringstream oss;
  oss << std::numeric_limits<int64_t>::max();
  return oss.str();
#endif
}

} // namespace primecount
