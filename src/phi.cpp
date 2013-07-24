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
/// This implementation is based on an algorithm from Tomas Oliveira e
/// Silva [1]. I have added a cache to my implementation in which
/// results of phi(x, a) are stored if x is small. This cache speeds
/// up the calculations by at least 3 orders of magnitude near 10^15.
/// For PHI_CACHE_LIMIT = 32767 the memory usage of the cache is < 200
/// megabytes per thread.
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
    int size = ps.countPrimes(0, PHI_CACHE_LIMIT);
    cache_.resize(size);
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
        // x2 = x / primes_[a2]
        int64_t x2 = FAST_DIV(x, primes_[a2]);

        if (x2 < PHI_CACHE_LIMIT && validIndexes(a2, x2) && 
            cache_[a2][x2] != 0)
        {
          // phi(x2, a2) is cached
          sum += cache_[a2][x2] * -SIGN;
        }
        else
        {
          // phi(x2, a2) is not cached, calculate recursively
          int64_t phiResult = phi<-SIGN>(x2, a2);
          // cache phi(x2, a2)
          if (x2 < PHI_CACHE_LIMIT)
          {
            if (!validIndexes(a2, x2))
              cache_[a2].resize(x2 + 1, 0);
            cache_[a2][x2] = static_cast<int16_t>(phiResult * -SIGN);
          }
          sum += phiResult;
        }
      }
    }
    return sum;
  }
private:
  /// First a primes needed to calculate phi(x, a)
  const std::vector<int32_t>& primes_;
  /// Cache of phi(x, a) results with x < PHI_CACHE_LIMIT
  std::vector<std::vector<int16_t> > cache_;

  bool validIndexes(int64_t a2, int64_t x2) const
  {
    return x2 < static_cast<int64_t>(cache_[a2].size());
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

} // namespace primecount
