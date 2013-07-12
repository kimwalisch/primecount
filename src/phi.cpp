///
/// @file  phi.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pi_bsearch.h"
#include "isqrt.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>
#include <limits>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.h"
#endif

/// Results of phi(x, a) are cached for x < PHI_CACHE_LIMIT
/// @pre PHI_CACHE_LIMIT <= 32767
///
#define PHI_CACHE_LIMIT 32767

/// Avoids slow 64-bit division if possible
#define FAST_DIV(x, y) ((x <= std::numeric_limits<uint32_t>::max()) \
    ? static_cast<uint32_t>(x) / (y) : (x) / (y))


namespace primecount {

/// This class is used to calculate phi(x, a) using the recursive
/// formula: phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1).
/// This implementation is based on an algorithm described in a
/// paper [1] by Tomas Oliveira e Silva.
/// On top of that algorithm I have implemented an efficient caching
/// solution that caches results of small (x, a) pairs in a two
/// dimensional vector. The memory usage is up to 200 megabytes (per
/// thread) for PHI_CACHE_LIMIT = 32767.
///
/// [1] Tomas Oliveira e Silva, "Computing pi(x): the combinatorial method",
///     Revista do DETUA, vol. 4, no. 6, pp. 759-768, March 2006
///
class PhiCache {
public:
  PhiCache(const std::vector<int32_t>& primes)
    : primes_(primes)
  {
    PrimeSieve ps;
    cache_.resize(ps.countPrimes(0, PHI_CACHE_LIMIT));
  }

  template<int64_t SIGN> int64_t phi(int64_t x, int64_t a)
  {
    int64_t sum = x * SIGN;
    if (a > 0)
    {
      int64_t iters = pi_bsearch(primes_.begin(), primes_.begin() + a, isqrt(x));
      sum += (a - iters) * -SIGN;

      for (int64_t a2 = 0; a2 < iters; a2++)
      {
        // next x = x / primes_[a2]
        int64_t x2 = FAST_DIV(x, primes_[a2]);
        int64_t phiValue;

        if (x2 < PHI_CACHE_LIMIT && x2 < getCacheSize(a2) && 
            cache_[a2][x2] != 0)
        {
          // phi(x2, a2) is cached
          phiValue = cache_[a2][x2] * -SIGN;
        }
        else
        {
          // phi(x2, a2) is not cached, calculate recursively
          phiValue = phi<-SIGN>(x2, a2);

          if (x2 < PHI_CACHE_LIMIT)
          {
            if (x2 >= getCacheSize(a2))
              cache_[a2].resize(x2 + 1, 0);

            // cache phi(x2, a2)
            cache_[a2][x2] = static_cast<int16_t>(phiValue * -SIGN);
          }
        }
        sum += phiValue;
      }
    }
    return sum;
  }
private:
  /// First a primes needed to calculate phi(x, a)
  const std::vector<int32_t>& primes_;

  /// Cache of phi(x, a) values for small values of x
  /// The memory usage of cache_ is:
  /// pi(PHI_CACHE_LIMIT) * PHI_CACHE_LIMIT * sizeof(int16_t)
  ///
  std::vector<std::vector<int16_t> > cache_;

  int64_t getCacheSize(int64_t a2) const
  {
    return static_cast<int64_t>(cache_[a2].size());
  }
};

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads)
{
  if (x < 1) return 0;
  if (a < 1) return x;

  std::vector<int32_t> primes;
  PrimeSieve ps;
  ps.generate_N_Primes(a, &primes);

  int iters = pi_bsearch(primes.begin(), primes.begin() + a, isqrt(x));
  PhiCache cache(primes);
  int64_t sum = x - a + iters;

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
      num_threads(threads) schedule(dynamic, 16)
#endif
  for (int i = 0; i < iters; i++)
    sum += cache.phi<-1>(x / primes[i], i);

  return sum;
}

/// 2nd partial sieve function.
/// P2(x, a) counts the numbers <= x that have exactly 2 prime
/// factors which are all greater than the a-th prime.
///
int64_t P2(int64_t x, int64_t a, int64_t b, int threads)
{
  int64_t sum = (b + a - 2) * (b - a + 1) / 2;
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
    sum -= pix;
  }

  return sum;
}

/// 3rd partial sieve function.
/// P3(x, a) counts the numbers <= x that have exactly 3 prime
/// factors which are all greater than the a-th prime.
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
      sum2 -= pi_bsearch(primes.begin(), primes.end(), x2 / primes[j - 1]) - (j - 1);

    sum += sum2;
  }

  return sum;
}

} // namespace primecount
