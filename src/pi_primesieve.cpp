#include <primecount.h>
#include <primesieve/soe/ParallelPrimeSieve.h>
#include <stdint.h>

namespace primecount {

int64_t pi_primesieve(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS)
    pps.setNumThreads( threads );

  int64_t sum = pps.countPrimes(0, x);

  return sum;
}

} // namespace primecount
