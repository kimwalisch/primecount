///
/// @file  phi.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pi_bsearch.h"
#include "imath.h"

#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>
#include <limits>
#include <cstddef>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.h"
#endif

/// Avoids slow 64-bit division if possible
#define FAST_DIV(x, y) ((x <= std::numeric_limits<uint32_t>::max()) \
    ? static_cast<uint32_t>(x) / (y) : (x) / (y))

namespace primecount {

/// This class is used to calculate phi(x, a) using the recursive
/// formula: phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1).
/// This implementation is based on an algorithm from Tomas Oliveira e
/// Silva [1]. I have added a cache to my implementation in which
/// results of phi(x, a) are stored if x < 2^16 and a < 500. This
/// cache speeds up the calculations by at least 3 orders of magnitude
/// near 10^15.
///
/// [1] Tomas Oliveira e Silva, "Computing pi(x): the combinatorial method",
///     Revista do DETUA, vol. 4, no. 6, pp. 759-768, March 2006
///
class PhiCache {
public:
  PhiCache(const std::vector<int32_t>& primes) :
    primes_(primes), bytes_(0)
  {
    std::size_t max_size = CACHE_A_LIMIT;
    cache_.resize(std::min(primes.size(), max_size));
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

        if (validIndexes(x2, a2) && cache_[a2][x2] != 0)
          sum += cache_[a2][x2] * -SIGN;
        else
        {
          // phi(x2, a2) is not cached
          int64_t phiResult = phi<-SIGN>(x2, a2);
          if (validIndexes(x2, a2))
            cache_[a2][x2] = static_cast<uint16_t>(phiResult * -SIGN);
          sum += phiResult;
        }
      }
    }
    return sum;
  }
private:
  enum {
    /// Keep the cache size below CACHE_BYTES_LIMIT per thread
    CACHE_BYTES_LIMIT = 16 << 20,
    CACHE_A_LIMIT = 500
  };
  /// Cache for phi(x, a) results
  std::vector<std::vector<uint16_t> > cache_;
  const std::vector<int32_t>& primes_;
  int64_t bytes_;

  bool validIndexes(int64_t x2, int64_t a2)
  {
    if (a2 >= CACHE_A_LIMIT || x2 >= std::numeric_limits<uint16_t>::max())
      return false;
    // resize cache if necessary
    if (x2 >= static_cast<int64_t>(cache_[a2].size()))
    {
      if (bytes_ > CACHE_BYTES_LIMIT)
        return false;
      bytes_ += (x2 + 1 - cache_[a2].size()) * 2;
      cache_[a2].resize(x2 + 1, 0);
    }
    return true;
  }
};

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
/// @pre prime[a] <= sqrt(x)
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
  for (int64_t a2 = 0; a2 < iters; a2++)
    sum += cache.phi<-1>(x / primes[a2], a2);

  return sum;
}

} // namespace primecount
