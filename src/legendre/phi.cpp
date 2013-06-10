#include "../utils/isqrt.h"
#include "../utils/Next_N_Primes_Vector.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
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

/// This method is an optimized version (binary search) of
/// int i = 0;
/// while (i < a && primes_[i] <= sqrt(x))
///   i++;
/// return i;
///
uint32_t findSqrtIndex(int64_t x, int64_t a, const std::vector<uint32_t>& primes)
{
  uint32_t root = isqrt(x);
  uint32_t index = static_cast<uint32_t>(
      std::upper_bound(primes.begin(), primes.begin() + a, root)
      - primes.begin());
  return index;
}

class Cache {
public:
  Cache(const std::vector<uint32_t>& primes)
    : primes_(primes)
  {
    PrimeSieve ps;
    phiCache_.resize(ps.countPrimes(0, PHI_CACHE_LIMIT));
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
      int64_t limit = findSqrtIndex(x, a, primes_);
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
  Next_N_Primes_Vector<uint32_t> primes;
  primes.generatePrimes(/* start = */ 0 , /* n = */ a);
  int limit = findSqrtIndex(x, a, primes);
  Cache cache(primes);
  int64_t sum = x - a + limit;

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
      num_threads(threads) schedule(static, 128)
#endif
  for (int i = 0; i < limit; i++)
    sum += cache.phi<-1>(x / primes[i], i);

  return sum;
}

} // namespace legendre
