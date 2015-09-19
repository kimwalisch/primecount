///
/// @file   primecount.cpp
/// @brief  primecount C++ API
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primecount.hpp>
#include <calculator.hpp>
#include <int128.hpp>
#include <pmath.hpp>

#include <algorithm>
#include <ctime>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace {

int threads_ = primecount::MAX_THREADS;

double alpha_ = -1;

}

namespace primecount {

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
/// @pre   x  <= get_max_x().
///
string pi(const string& x)
{
  return pi(x, threads_);
}

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like "10^12".
/// @pre   x  <= get_max_x().
///
string pi(const string& x, int threads)
{
  maxint_t n = to_maxint(x);
  maxint_t pin = pi(n, threads);
  ostringstream oss;
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
  return pi_deleglise_rivat_parallel2(x, threads);
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
  if (x <= numeric_limits<int64_t>::max())
    return pi_deleglise_rivat((int64_t) x, threads);

  return pi_deleglise_rivat_parallel3(x, threads);
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
/// pi(string x). The return type is a string as max may be a
/// 128-bit integer which is not supported by all compilers.
///
string get_max_x(double alpha)
{
  ostringstream oss;

#ifdef HAVE_INT128_T
  // primecount is limited by:
  // z < 2^62, with z = x^(2/3) / alpha
  // x^(2/3) / alpha < 2^62
  // x < (2^62 * alpha)^(3/2)

  // safety buffer: use 61 instead of 62 
  double max_x = pow(pow(2.0, 61.0) * alpha, 3.0 / 2.0);
  oss << (int128_t) max_x; 
#else
  unused_param(alpha); 
  oss << numeric_limits<int64_t>::max();
#endif

  return oss.str();
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
  threads = (int) min((int64_t) threads, sieve_limit / thread_threshold);
  threads = max(1, threads);
  return threads;
}

void set_alpha(double alpha)
{
  alpha_ = alpha;
}

double get_alpha()
{
  return alpha_;
}

double get_alpha(maxint_t x, int64_t y)
{
  // y = x13 * alpha, thus alpha = y / x13
  double x13 = (double) iroot<3>(x);
  return (double) y / x13;
}

/// Calculate the Lagarias-Miller-Odlyzko alpha tuning factor.
/// alpha = a log(x)^2 + b log(x) + c
/// a, b and c are constants that should be determined empirically.
/// @see ../doc/alpha-factor-tuning.pdf
///
double get_alpha(maxint_t x, double a, double b, double c)
{
  double alpha = get_alpha();

  if (alpha < 1)
  {
    double logx = log((double) x);
    alpha = a * pow(logx, 2) + b * logx + c;
  }

  return in_between(1, alpha, iroot<6>(x));
}

/// Calculate the Deleglise-Rivat alpha tuning factor.
/// alpha = a log(x)^3 + b log(x)^2 + c log(x) + d
/// a, b, c and d are constants that should be determined empirically.
/// @see ../doc/alpha-factor-tuning.pdf
///
double get_alpha(maxint_t x, double a, double b, double c, double d)
{
  double alpha = get_alpha();

  if (alpha < 1)
  {
    double logx = log((double) x);
    alpha = a * pow(logx, 3) + b * pow(logx, 2) + c * logx + d;
  }

  return in_between(1, alpha, iroot<6>(x));
}

void set_num_threads(int threads)
{
  threads_ = validate_threads(threads);
}

int get_num_threads()
{
  return validate_threads(threads_);
}

maxint_t to_maxint(const string& expr)
{
  maxint_t n = calculator::eval<maxint_t>(expr);
  return n;
}

} // namespace primecount
