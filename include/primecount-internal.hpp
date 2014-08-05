///
/// @file   primecount-internal.hpp
/// @brief  primecount internal function definitions.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_INTERNAL_HPP
#define PRIMECOUNT_INTERNAL_HPP

#include <ptypes.hpp>
#include <aligned_vector.hpp>

#include <stdint.h>
#include <string>
#include <vector>

namespace primecount {

enum {
  /// Uses all CPU cores.
  MAX_THREADS = -1
};

class PhiCache;

int64_t pi(int64_t x, int threads);

#ifdef HAVE_INT128_T

int128_t pi(int128_t x);

int128_t pi(int128_t x, int threads);

#endif

/// Alias for the fastest prime counting function in primecount.
/// @param x  integer or arithmetic expression like 10^12.
/// @pre   x  <= primecount::max().
///
std::string pi(const std::string& x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat(int64_t x, int threads);

#ifdef HAVE_INT128_T

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat(int128_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat(int128_t x, int threads);

#endif

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat1(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat2(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat3(int64_t x);

#ifdef HAVE_INT128_T

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat4(int128_t x);

#endif

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel3(int64_t x, int threads);

#ifdef HAVE_INT128_T

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat_parallel4(int128_t x, int threads);

#endif

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^(1/2)) space.
///
int64_t pi_legendre(int64_t x, int threads);

/// Calculate the number of primes below x using Lehmer's formula.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2)) space.
///
int64_t pi_lehmer(int64_t x, int threads);

/// Calculate the number of primes below x using Lehmer's formula.
/// This version uses a different P2(x, y) implementation,
/// it runs slower than pi_lehmer(x) on most systems.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2) / log x) space.
///
int64_t pi_lehmer2(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses the
/// recursive phi formula with caching to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
int64_t pi_lmo1(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(2/3) / (log x)^2) space.
///
int64_t pi_lmo2(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo3(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes and a special data structure
/// for counting to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo4(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo5(int64_t x);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo_parallel1(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo_parallel2(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo_parallel3(int64_t x, int threads);

/// Calculate the number of primes below x using Meissel's formula.
/// Run time: O(x/(log x)^3) operations, O(x^(1/2) / log x) space.
///
int64_t pi_meissel(int64_t x, int threads);

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x, int threads);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, PhiCache* phiCache);

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

/// Calculate the nth prime using a combination of an efficient prime
/// counting function implementation and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/2)) space.
///
int64_t nth_prime(int64_t n, int threads);

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O((x / primes[a])^(1/2)).
///
int64_t P2(int64_t x, int64_t y, int threads);

#ifdef HAVE_INT128_T

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O((x / primes[a])^(1/2)).
///
int128_t P2(int128_t x, int64_t y, int threads);

#endif

/// P2_lehmer(x, a) counts the numbers <= x that have exactly 2
/// prime factors each exceeding the a-th prime.
/// Space complexity: O(pi(x^(1/2))).
///
int64_t P2_lehmer(int64_t x, int64_t a, int threads);

/// P3(x, a) counts the numbers <= x that have exactly 3
/// prime factorseach exceeding the a-th prime.
/// Space complexity: O(pi(x^(1/2))).
///
int64_t P3(int64_t x, int64_t a, int threads);

/// Used to balance the load in the computation of the special
/// leaves in the LMO and Deleglise-Rivat algorithms.
///
void balance_S2_load(double x,
                     double threads,
                     double* old_rsd,
                     aligned_vector<double>& timings,
                     int64_t* segment_size,
                     int64_t* segments_per_thread,
                     int64_t min_segment_size,
                     int64_t max_segment_size);

/// Convert a string into an integer of type maxint_t.
maxint_t to_maxint(const std::string& expr);

/// Get the wall time in seconds.
double get_wtime();

int validate_threads(int threads);

int validate_threads(int threads, int64_t sieve_limit, int64_t thread_threshold = 100000);

} // namespace primecount

#endif
