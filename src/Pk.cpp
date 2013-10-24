///
/// @file  Pk.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pi_bsearch.h"
#include "imath.h"

#include <primesieve.hpp>
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
/// Space complexity: O(x^0.5).
///
int64_t P2(int64_t x, int64_t a, int threads)
{
  std::vector<int32_t> primes;
  std::vector<int64_t> counts;
  primesieve::generate_primes(isqrt(x), &primes);
  counts.resize(primes.size());

  int64_t b = pi_bsearch(primes.begin(), primes.end(), isqrt(x));
  int64_t sum = 0;
  int64_t pix = 0;

  // This uses a clever trick, instead of calculating
  // pi(x / primes[i]) for a < i <= b it only counts the primes
  // between adjacent values [x/primes[i], x/primes[i-1]].
  // When finished pi(x / primes[i]) can quickly be calculated
  // by backwards summing up the counts.
  //
#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = b; i > a; i--)
  {
    int64_t prev = (i != b) ? x / primes[i] + 1 : 0;
    int64_t xi = x / primes[i - 1];
    counts[i - 1] = primesieve::count_primes(prev, xi);
  }

  for (int64_t i = b; i > a; i--)
  {
    pix += counts[i - 1];
    sum += pix;
    sum -= i - 1;
  }

  return sum;
}

/// 3rd partial sieve function.
/// P3(x, a) counts the numbers <= x that have exactly 3 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O(x^0.5).
///
int64_t P3(int64_t x, int64_t a, int threads)
{
  std::vector<int32_t> primes;
  primesieve::generate_primes(isqrt(x), &primes);

  int64_t c = pi_bsearch(primes.begin(), primes.end(), iroot<3>(x));
  int64_t sum = 0;

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for reduction(+: sum) schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = a + 1; i <= c; i++)
  {
    int64_t xi = x / primes[i - 1];
    int64_t bi = pi_bsearch(primes.begin(), primes.end(), isqrt(xi));

    for (int64_t j = i; j <= bi; j++)
    {
      int64_t xij = xi / primes[j - 1];
      sum += pi_bsearch(primes.begin(), primes.end(), xij) - (j - 1);
    }
  }

  return sum;
}

} // namespace primecount
