///
/// @file  Pk.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pi_bsearch.h"
#include "isqrt.h"

#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.h"
#endif

namespace primecount {

/// 2nd partial sieve function.
/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime.
///
int64_t P2(int64_t x, int64_t a, int64_t b, int threads)
{
  int64_t sum = -((b + a - 2) * (b - a + 1) / 2);
  int64_t pix = 0;
  std::vector<int32_t> primes;
  std::vector<int64_t> counts;

  PrimeSieve ps;
  ps.generate_N_Primes(b, &primes);
  counts.resize(primes.size());

  // This uses a clever trick, instead of calculating
  // pi(x / primes[i]) for a < i <= b it only counts the primes
  // between adjacent values [x / primes[i], x / primes[i - 1]].
  // When finished pi(x / primes[i]) can quickly be calculated
  // by backwards summing up the counts.
  //
#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for private(ps) schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = b; i > a; i--)
  {
    int64_t x2 = (i != b) ? x / primes[i] + 1 : 0;
    int64_t x3 = x / primes[i - 1];
    counts[i - 1] = ps.countPrimes(x2, x3);
  }

  for (int64_t i = b; i > a; i--)
  {
    pix += counts[i - 1];
    sum += pix;
  }

  return sum;
}

/// 3rd partial sieve function.
/// P3(x, a) counts the numbers <= x that have exactly 3 prime
/// factors each exceeding the a-th prime.
///
int64_t P3(int64_t x, int64_t a, int64_t b, int64_t c, int threads)
{
  std::vector<int32_t> primes;
  PrimeSieve ps;
  ps.generate_N_Primes(b, &primes);
  int64_t sum = 0;

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for reduction(+: sum) schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = a + 1; i <= c; i++)
  {
    int64_t x2 = x / primes[i - 1];
    int64_t bi = pi_bsearch(primes.begin(), primes.end(), isqrt(x2));
    int64_t sum2 = 0;

    for (int64_t j = i; j <= bi; j++)
      sum2 += pi_bsearch(primes.begin(), primes.end(), x2 / primes[j - 1]) - (j - 1);

    sum += sum2;
  }

  return sum;
}

} // namespace primecount
