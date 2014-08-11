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
#include <pmath.hpp>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace {

// By default primecount uses all CPU cores.
int threads_ = primecount::MAX_THREADS;

// Print status information to std::cout.
bool print_status_ = false;

}

namespace primecount {

void set_print_status(bool print_status)
{
  print_status_ = print_status;
}

bool print_status()
{
  return print_status_;
}

void print_seconds(double seconds)
{
  std::cout << "Seconds: " << std::fixed << std::setprecision(3) << seconds << std::endl;
}

void print_megabytes(std::size_t bytes)
{
  double megabytes = bytes / (double) (1 << 20);
  std::cout << "memory usage = " << std::fixed << std::setprecision(3) << megabytes << " megabytes" << std::endl;
}

/// Get the wall time in seconds.
double get_wtime()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
#endif
}

int validate_threads(int threads)
{
#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  return in_between(1, threads, omp_get_max_threads());
#else
  threads = 1;
  return threads; 
#endif
}

int validate_threads(int threads, int64_t sieve_limit, int64_t thread_threshold)
{
  threads = validate_threads(threads);
  threads = (int) std::min((int64_t) threads, sieve_limit / thread_threshold);
  threads = std::max(1, threads);
  return threads;
}

void set_num_threads(int threads)
{
  threads_ = validate_threads(threads);
}

int get_num_threads()
{
  return validate_threads(threads_);
}

maxint_t to_maxint(const std::string& expr)
{
  maxint_t n = calculator::eval<maxint_t>(expr);
  return n;
}

int64_t pi(int64_t x)
{
  return pi(x, threads_);
}

int64_t pi(int64_t x, int threads)
{
  return pi_deleglise_rivat(x, threads);
}

#ifdef HAVE_INT128_T

int128_t pi(int128_t x)
{
  return pi(x, threads_);
}

int128_t pi(int128_t x, int threads)
{
  return pi_deleglise_rivat(x, threads);
}

#endif

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like "10^12".
/// @pre   x  <= primecount::max().
///
std::string pi(const std::string& x)
{
  return pi(x, threads_);
}

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like "10^12".
/// @pre   x  <= primecount::max().
///
std::string pi(const std::string& x, int threads)
{
  maxint_t n = to_maxint(x);
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
  if (threads <= 1)
    return pi_deleglise_rivat3(x);

  return pi_deleglise_rivat_parallel3(x, threads);
}

#ifdef HAVE_INT128_T

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat(int128_t x)
{
  return pi_deleglise_rivat(x, threads_);
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat(int128_t x, int threads)
{
  // use 64-bit if possible
  if (x <= std::numeric_limits<int64_t>::max())
    return pi_deleglise_rivat((int64_t) x, threads);

  if (threads <= 1)
    return pi_deleglise_rivat4(x);
  else
    return pi_deleglise_rivat_parallel4(x, threads);
}

#endif

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
  if (threads <= 1)
    return pi_lmo5(x);

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
  // max = 10^27 if int128_t is supported
  return std::string("1") + std::string(27, '0');
#else
  std::ostringstream oss;
  oss << std::numeric_limits<int64_t>::max();
  return oss.str();
#endif
}

} // namespace primecount
