///
/// @file  primecount.hpp
/// @brief primecount C++ API. primecount is a C/C++ library for
///        counting the number of primes <= x using highly
///        optimized implementations of the combinatorial type
///        prime counting function algorithms.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License.
///

#ifndef PRIMECOUNT_HPP
#define PRIMECOUNT_HPP

#include <stdexcept>
#include <string>
#include <stdint.h>

#define PRIMECOUNT_VERSION "7.20"
#define PRIMECOUNT_VERSION_MAJOR 7
#define PRIMECOUNT_VERSION_MINOR 20

namespace primecount {

class primecount_error : public std::runtime_error
{
public:
  primecount_error(const std::string& msg)
    : std::runtime_error(msg)
  { }
};

/// pc_int128_t is a portable int128_t type used by primecount's C++ API.
///
/// How to convert a pc_int128_t to an int128_t:
/// int128_t n = x.lo | (int128_t(x.hi) << 64);
///
/// How to convert an int128_t to a pc_int128_t:
/// pc_int128_t x;
/// x.lo = (uint64_t) n;
/// x.hi = (int64_t) (n >> 64);
///
struct pc_int128_t
{
  uint64_t lo;
  int64_t hi;
};

/// Count the number of primes <= x using Xavier Gourdon's
/// algorithm. Uses all CPU cores by default.
/// Throws a primecount_error if an error occurs.
///
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi(int64_t x);

/// 128-bit prime counting function.
/// Count the number of primes <= x using Xavier Gourdon's
/// algorithm. Uses all CPU cores by default.
/// Throws a primecount_error if an error occurs.
/// 
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
pc_int128_t pi(pc_int128_t x);

/// 128-bit prime counting function.
/// Count the number of primes <= x using Xavier Gourdon's
/// algorithm. Uses all CPU cores by default.
///
/// @param x Null-terminated string integer e.g. "12345".
///          Note that x must be <= get_max_x() which is 10^31 on
///          64-bit systems and 2^63-1 on 32-bit systems.
/// Throws a primecount_error if an error occurs.
///
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
std::string pi(const std::string& x);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
/// Throws a primecount_error if an error occurs.
///
int64_t phi(int64_t x, int64_t a);

/// Find the nth prime using a combination of the prime counting
/// function and the sieve of Eratosthenes.
/// @pre n <= 216289611853439384
/// Throws a primecount_error if an error occurs.
///
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t nth_prime(int64_t n);

/// 128-bit nth prime function.
/// Find the nth prime using a combination of the prime counting
/// function and the sieve of Eratosthenes.
///
/// @pre n <= 10^29 on 64-bit systems and
///      n <= 216289611853439384 on 32-bit systems.
/// Throws a primecount_error if an error occurs.
///
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
pc_int128_t nth_prime(pc_int128_t n);

/// Largest number supported by pi(const std::string& x).
/// @return 64-bit CPUs: 10^31,
///         32-bit CPUs: 2^63-1.
/// The return type is a string as get_max_x() may be a 128-bit
/// integer which is not supported by some compilers.
///
std::string get_max_x();

/// Get the currently set number of threads
int get_num_threads();

/// Set the number of threads
void set_num_threads(int num_threads);

/// Recompute pi(x) with alternative alpha tuning factor(s) to
/// verify the first result. This redundancy helps guard 
/// against potential bugs in primecount: if an error exists,
/// it is highly unlikely that both pi(x) computations would
/// produce the same (incorrect) result.
///
void set_verify_computation(bool enable);

/// Get the primecount version number, in the form “i.j”
std::string primecount_version();

} // namespace

#endif
