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

#include <primecount.hpp>
#include <stdint.h>
#include <vector>

namespace primecount {

class PhiCache;

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat1(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat2(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using Lehmer's formula.
/// This version uses a different P2(x, y) implementation,
/// it runs slower than pi_lehmer(x) on most systems.
/// Run time: O(x/(log x)^4) operations, O(x^(1/2) / log x) space.
///
int64_t pi_lehmer2(int64_t x, int threads = MAX_THREADS);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses the
/// recursive phi formula with caching to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo1(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(2/3) / log log x) space.
///
int64_t pi_lmo2(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo3(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm. This implementation uses
/// the segmented sieve of Eratosthenes and a special data structure
/// for counting to calculate S2(x).
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo4(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo5(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo6(int64_t x);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel1(int64_t x, int threads = MAX_THREADS);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel2(int64_t x, int threads = MAX_THREADS);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel3(int64_t x, int threads = MAX_THREADS);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel4(int64_t x, int threads = MAX_THREADS);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, PhiCache* phiCache);

/// Calculates the number of 1 bits inside an array.
int64_t popcount(const uint64_t* bits, int64_t start, int64_t stop);

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O((x / primes[a])^(1/2)).
///
int64_t P2(int64_t x, int64_t y, int threads = MAX_THREADS);

/// P2_lehmer(x, a) counts the numbers <= x that have exactly 2
/// prime factors each exceeding the a-th prime.
/// Space complexity: O(pi(x^(1/2))).
///
int64_t P2_lehmer(int64_t x, int64_t a, int threads = MAX_THREADS);

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
