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

int64_t P2(int64_t x, int64_t a, int64_t b, int64_t pb, int threads /* = MAX_THREADS */)
{
  PrimeSieveVector<uint32_t> primes;
  PrimeSieve ps;
  ps.generatePrimes(0, pb, &primes);

  int64_t sum = (b + a - 2) * (b - a + 1) / 2;
  int64_t pix = 0;
  int64_t old = 0;

  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS)
    pps.setNumThreads( threads );

  for (int i = static_cast<int>(primes.size()); i > a; i--)
  {
    int64_t x2 = x / primes[i - 1];
    pix += pps.countPrimes(old + 1, x2);
    old = x2;
    sum -= pix;
  }

  return sum;
}

} // namespace primecount
