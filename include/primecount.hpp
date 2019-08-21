///
/// @file  primecount.hpp
/// @brief primecount C++ API
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License.
///

#ifndef PRIMECOUNT_HPP
#define PRIMECOUNT_HPP

#include <stdexcept>
#include <string>
#include <stdint.h>

#define PRIMECOUNT_VERSION "5.0"
#define PRIMECOUNT_VERSION_MAJOR 5
#define PRIMECOUNT_VERSION_MINOR 0

namespace primecount {

class primecount_error : public std::runtime_error
{
public:
  primecount_error(const std::string& msg)
    : std::runtime_error(msg)
  { }
};

/// Count the number of primes <= x using Xavier Gourdon's
/// algorithm. Uses all CPU cores by default.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi(int64_t x);

/// Count the number of primes <= x using Xavier Gourdon's
/// algorithm. Uses all CPU cores by default.
/// @param x Number or arithmetic expression e.g. "1000", "10^22"
/// @pre x <= 10^31 on 64-bit systems, x < 2^63 on 32-bit systems.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
std::string pi(const std::string& x);

/// Find the nth prime using a combination of the prime
/// counting function and the sieve of Eratosthenes.
/// @pre n <= 216289611853439384
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/2))
///
int64_t nth_prime(int64_t n);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a);

/// Set the number of threads
void set_num_threads(int num_threads);

/// Get the currently set number of threads
int get_num_threads();

/// Largest number supported by pi(const std::string& x).
/// @param alpha Tuning factor
/// @return 64-bit CPUs: max >= 10^27,
///         32-bit CPUs: 2^63-1
/// The return type is a string as get_max_x() may be a
/// 128-bit integer which is not supported by some
/// compilers.
///
std::string get_max_x(double alpha = 1.0);

/// Get the primecount version number, in the form “i.j”
std::string primecount_version();

} // namespace

#endif
