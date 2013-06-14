#include "utils/PrimeSieveVector.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <primesieve/soe/ParallelPrimeSieve.h>
#include <vector>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

int64_t P2(int64_t x, int64_t a, int64_t b, int64_t pb_max, int threads /* = MAX_THREADS */)
{
  PrimeSieveVector<uint32_t> primes;
  PrimeSieve ps;
  ps.generatePrimes(0, pb_max, &primes);

  int64_t sum = (b + a - 2) * (b - a + 1) / 2;
  int64_t pix = 0;
  int64_t old = 0;

  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS)
    pps.setNumThreads( threads );

  int i = static_cast<int>(primes.size());

  // optimized for large intervals
  for (; i > a; i--)
  {
    int64_t x2 = x / primes[i - 1];
    if (x2 - old < 100000000)
      break;
    pix += pps.countPrimes(old + 1, x2);
    old = x2;
    sum -= pix;
  }

  // optimized for small intervals
#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for private(ps) firstprivate(pix, old) reduction(+: sum) \
      num_threads(threads) schedule(dynamic)
#endif
  for (int j = i; j > a; j--)
  {
    int64_t x2 = x / primes[j - 1];
    pix += ps.countPrimes(old + 1, x2);
    old = x2;
    sum -= pix;
  }

  return sum;
}

} // namespace primecount
