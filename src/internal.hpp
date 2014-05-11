///
/// @file   internal.hpp
/// @brief  primecount internal function definitions.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_INTERNAL_HPP
#define PRIMECOUNT_INTERNAL_HPP

#include <primecount.hpp>
#include <stdint.h>
#include <vector>

namespace primecount {

class PhiCache;

/// Calculate the number of primes below x using Lehmer's formula.
/// This version is for testing only, it uses a different P2(x, a)
/// implementation than pi_lehmer(x).
///
int64_t pi_lehmer2(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses the
/// recursive phi formula with caching to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo1(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(2/3) / log log x) space.
///
int64_t pi_lmo2(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo3(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes and a special data structure
/// for counting to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo4(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo5(int64_t x, int threads = MAX_THREADS);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, PhiCache* phiCache);

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O((x / primes[a])^(1/2)).
///
int64_t P2(int64_t x, int64_t a);

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O(pi(x^(1/2))).
///
int64_t P2(int64_t x, int64_t a, int threads);

/// P3(x, a) counts the numbers <= x that have exactly 3
/// prime factorseach exceeding the a-th prime.
/// Space complexity: O(pi(x^(1/2))).
///
int64_t P3(int64_t x, int64_t a, int threads = MAX_THREADS);

/// Calculate the contribution of the ordinary leaves in the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(y) operations, O(y) space.
///
int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           std::vector<int32_t>& primes,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu);

} // namespace primecount

#endif
