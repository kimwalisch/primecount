///
/// @file  primecount.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License.
///

#ifndef PRIMECOUNT_HPP
#define PRIMECOUNT_HPP

#include <stdexcept>
#include <string>
#include <stdint.h>

#define PRIMECOUNT_VERSION "3.5"
#define PRIMECOUNT_VERSION_MAJOR 3
#define PRIMECOUNT_VERSION_MINOR 5

namespace primecount {

class primecount_error : public std::runtime_error
{
public:
  primecount_error(const std::string& msg)
    : std::runtime_error(msg)
  { }
};

/// Alias for the fastest prime counting function in primecount.
int64_t pi(int64_t x);

/// 128-bit prime counting function.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
/// @param expr  Integer arithmetic expression e.g. "1000", "10^22"
/// @pre   expr  <= get_max_x()
///
std::string pi(const std::string& x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
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
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
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

/// Calculate the nth prime using a combination of the counting
/// function and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/2)) space.
/// @pre n <= 216289611853439384
///
int64_t nth_prime(int64_t n);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a);

/// Calculate the offset logarithmic integral which is a very accurate
/// approximation of the number of primes below x.
/// @post Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
int64_t Li(int64_t x);

/// Calculate the inverse logarithmic integral Li^-1(x) which is
/// a very accurate approximation of the nth prime.
/// @post Li_inverse(x) < nth_prime(x) for 7 <= x <= ~ 10^316
///
int64_t Li_inverse(int64_t x);

/// Enable/disable printing status information during computation.
void set_print_status(bool print_status);

/// Set the number of threads.
void set_num_threads(int num_threads);

/// Get the currently set number of threads.
int get_num_threads();

/// Largest integer supported by pi(const std::string& x).
/// The return type is a string as max may be a 128-bit integer
/// which is not supported by some compilers.
/// @param alpha Tuning factor used in LMO type algorithms.
/// @return for 32-bit CPUs: 2^63-1, 
///         for 64-bit CPUs: max >= 10^27
///
std::string get_max_x(double alpha = 1.0);

/// Get the primecount version number, in the form “i.j”.
std::string primecount_version();

} // namespace

#endif
