#include "phi.h"
#include "../utils/isqrt.h"
#include "../utils/PrimeSieveVector.h"

#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>
#include <limits>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// Results of phi(x, a) are cached for x < PHI_CACHE_LIMIT
/// @pre PHI_CACHE_LIMIT < 65535
///
#define PHI_CACHE_LIMIT 32767

/// Avoids slow 64-bit division if possible
#define FAST_DIV(x, y) ((x <= std::numeric_limits<uint32_t>::max()) \
  ? static_cast<uint32_t>(x) / (y) : (x) / (y))


namespace legendre {

class Cache {
public:
  Cache(const std::vector<uint32_t>& primes)
    : primes_(primes)
  {
    PrimeSieve ps;
    phiCache_.resize(ps.countPrimes(2, PHI_CACHE_LIMIT));
  }

  /// Calculate phi(x, a) using the recursive formula
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1).
  /// This implementation caches results of phi(x, a) for
  /// small values of x to speed up the calculations.
  ///
  template<int64_t SIGN> int64_t phi(int64_t x, int64_t a)
  {
    int64_t sum = x * SIGN;
    if (a > 0)
    {
      int64_t limit = findSqrtIndex(primes_, x, a);
      sum += (a - limit) * -SIGN;

      for (int64_t a2 = 0; a2 < limit; a2++)
      {
        // next x = x / primes_[a2]
        int64_t x2 = FAST_DIV(x, primes_[a2]);
        int64_t phiValue;

        if (x2 < PHI_CACHE_LIMIT && x2 < getCacheSize(a2) && 
            phiCache_[a2][x2] != 0)
        {
          // phi(x2, a2) is cached
          phiValue = phiCache_[a2][x2] * -SIGN;
        }
        else
        {
          // phi(x2, a2) is not cached, calculate recursively
          phiValue = phi<-SIGN>(x2, a2);

          if (x2 < PHI_CACHE_LIMIT)
          {
            if (x2 >= getCacheSize(a2))
              phiCache_[a2].resize(x2 + 1, 0);

            // cache phi(x2, a2)
            phiCache_[a2][x2] = static_cast<uint16_t>(phiValue * -SIGN);
          }
        }
        sum += phiValue;
      }
    }
    return sum;
  }
private:
  /// First a primes needed to calculate phi(x, a)
  const std::vector<uint32_t>& primes_;

  /// Cache of phi(x, a) values for small values of x
  /// The memory usage of phiCache_ is:
  /// pi(PHI_CACHE_LIMIT) * PHI_CACHE_LIMIT * sizeof(uint16_t)
  ///
  std::vector<std::vector<uint16_t> > phiCache_;

  int64_t getCacheSize(int64_t a2) const
  {
    return static_cast<int64_t>(phiCache_[a2].size());
  }
};

int64_t phi(int64_t x, int64_t a, int threads /* = MAX_THREADS */)
{
  int64_t sum = x;
  if (a > 0)
  {
    // store primes up to sqrt(x)
    PrimeSieveVector<uint32_t> primes;
    PrimeSieve ps;
    ps.generatePrimes(0, isqrt(x), &primes);
    primes.push_back(0);

    int limit = findSqrtIndex(primes, x, a);
    sum += (a - limit) * -1;
    Cache cache(primes);
#ifdef _OPENMP
    if (threads == MAX_THREADS)
      threads = omp_get_max_threads();
    #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
        num_threads(threads) schedule(static, 128)
#endif
    for (int a = 0; a < limit; a++)
      sum += cache.phi<-1>(x / primes[a], a);
  }
  return sum;
}

} // namespace legendre
