///
/// @file  pi_lmo1.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        Usually in the Lagarias-Miller-Odlyzko algorithm phi(x, a)
///        is calculated using a prime sieve but this simple
///        implementation calculates phi(x, a) using the recursive
///        formula with caching.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.hpp"
#include "PhiCache.hpp"
#include "pmath.hpp"

#include <primecount.hpp>
#include <primesieve.hpp>
#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
  #include "get_omp_threads.hpp"
#endif

using namespace std;

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo1(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x); 
  int64_t pi_y = pi_lehmer(y);
  int64_t c = min(PhiTiny::MAX_A, pi_y);
  int64_t S1 = 0;
  int64_t S2 = 0;

  vector<int32_t> lpf = make_least_prime_factor(y);
  vector<int32_t> mu = make_moebius(y);
  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(y, &primes);

  // Calculate the contribution of the ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += mu[n] * phi(x / n, c);

  const PhiTiny phiTiny;
  PhiCache cache(primes, phiTiny);

  // Calculate the contribution of the special leaves
#ifdef _OPENMP
  #pragma omp parallel for firstprivate(cache) schedule(dynamic) reduction(-: S2) \
      num_threads(get_omp_threads(threads)) 
#endif
  for (int64_t b = c + 1; b < pi_y; b++)
    for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= mu[m] * phi(x / (primes[b] * m), b - 1, &cache);

  int64_t phi = S1 + S2;
  int64_t sum = phi + pi_y - 1 - P2(x, pi_y);

  return sum;
}

} // namespace primecount
