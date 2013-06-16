#include "PrimeSieveVector.h"
#include "pi_bsearch.h"
#include "isqrt.h"

#include <primecount.h>
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


namespace primecount {

class PhiCache {
public:
  PhiCache(const std::vector<uint32_t>& primes)
    : primes_(primes)
  {
    cache_.resize(pi_primesieve(PHI_CACHE_LIMIT));
  }

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1).
  /// This implementation caches results of phi(x, a) for
  /// small values of x to speed up the calculations.
  ///
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
            cache_[a2][x2] = static_cast<uint16_t>(phiValue * -SIGN);
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
  /// The memory usage of cache_ is:
  /// pi(PHI_CACHE_LIMIT) * PHI_CACHE_LIMIT * sizeof(uint16_t)
  ///
  std::vector<std::vector<uint16_t> > cache_;

  int64_t getCacheSize(int64_t a2) const
  {
    return static_cast<int64_t>(cache_[a2].size());
  }
};

int64_t phi(int64_t x, int64_t a, int threads /* = MAX_THREADS */)
{
  PrimeSieveVector<uint32_t> primes;
  primes.generate_N_Primes(/* start = */ 0 , /* n = */ a);

  int iters = pi_bsearch(primes.begin(), primes.begin() + a, isqrt(x));
  PhiCache cache(primes);
  int64_t sum = x - a + iters;

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
      num_threads(threads) schedule(dynamic, 16)
#endif
  for (int i = 0; i < iters; i++)
    sum += cache.phi<-1>(x / primes[i], i);

  return sum;
}

} // namespace primecount
