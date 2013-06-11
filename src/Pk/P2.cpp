#include "../utils/isqrt.h"
#include "../utils/PrimeSieveVector.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

int64_t P2(int64_t x, int64_t a, int64_t pa, int64_t pb, int threads)
{
  PrimeSieveVector<uint32_t> primes;
  PrimeSieve ps;
  if (pa < pb)
    ps.generatePrimes(pa + 1, pb, &primes);

  int64_t b = a + static_cast<int64_t>( primes.size() );
  int64_t sum = (b + a - 2) * (b - a + 1) / 2;
  int64_t pix = 0;
  int64_t old = 0;

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for private(ps) firstprivate(pix, old) reduction(+: sum) \
      num_threads(threads) schedule(dynamic)
#endif
  for (int i = static_cast<int>(primes.size()) - 1; i >= 0; i--)
  {
    int64_t x2 = x / primes[i];
    if (old < x2)
      pix += ps.countPrimes(old + 1, x2);
    old = x2;
    sum -= pix;
  }

  return sum;
}

} // namespace primecount
