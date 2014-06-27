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

#include <stdint.h>
#include <vector>

namespace primecount {

enum {
  /// Uses all CPU cores.
  MAX_THREADS = -1
};

class PhiCache;
class FactorTable;

/// This is an alias for the fastest prime counting implementation
/// within primecount.
///
int64_t pi(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
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
int64_t pi_deleglise_rivat3(int64_t x);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads);

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_deleglise_rivat_parallel3(int64_t x, int threads);

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
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo(int64_t x, int threads);

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
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo5(int64_t x);

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo6(int64_t x);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel1(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel2(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel3(int64_t x, int threads);

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime counting algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel4(int64_t x, int threads);

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
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t nth_prime(int64_t n, int threads);

/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O((x / primes[a])^(1/2)).
///
int64_t P2(int64_t x, int64_t y, int threads);

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

/// Calculate the contribution of the ordinary leaves in the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(y) operations, O(y) space.
///
int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           std::vector<int32_t>& primes,
           FactorTable& Factors);

} // namespace primecount

#endif
