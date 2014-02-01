///
/// @file  pi_lmo_simple.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        Usually in the Lagarias-Miller-Odlyzko algorithm phi(x, a)
///        is calculated using a prime sieve but this simple
///        implementation calculates phi(x, a) using the recursive
///        formula with caching.
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.hpp"
#include "PhiCache.hpp"
#include "imath.hpp"
#include "Pk.hpp"

#include <primecount.hpp>
#include <primesieve.hpp>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.hpp"
#endif

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
/// @note O(x^0.5) space is due to parallel P2(x, a).
///
int64_t pi_lmo_simple(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t x13 = iroot<3>(x); 
  int64_t a = pi_lehmer(x13);
  int64_t a_minus_1 = a - 1;
  int64_t c = (a < 6) ? a : 6;
  int64_t S1 = 0;
  int64_t S2 = 0;

  std::vector<int32_t> primes;
  std::vector<int32_t> lpf = make_least_prime_factor(x13);
  std::vector<int32_t> mu = make_moebius(x13);

  primes.push_back(0);
  primesieve::generate_n_primes(a, &primes);

  const PhiTiny phiTiny;
  PhiCache cache(primes, phiTiny);

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for firstprivate(cache) reduction(+: S1) \
      num_threads(threads)
#endif
  for (int64_t n = 1; n <= x13; n++)
    if (lpf[n] > primes[c])
      S1 += mu[n] * phi(x / n, c, &cache);

#ifdef _OPENMP
  #pragma omp parallel for firstprivate(cache) reduction(-: S2) \
      num_threads(threads) schedule(dynamic)
#endif
  for (int64_t b = c; b < a_minus_1; b++)
    for (int64_t m = (x13 / primes[b + 1]) + 1; m <= x13; m++)
      if (lpf[m] > primes[b + 1])
        S2 -= mu[m] * phi(x / (m * primes[b + 1]), b, &cache);

  int64_t phi = S1 + S2;
  int64_t sum = phi + a - 1 - P2(x, a, threads);

  return sum;
}

} // namespace primecount
