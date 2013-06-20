#include "pi_bsearch.h"
#include "isqrt.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

int64_t P3(int64_t x, int64_t a, int64_t c, int64_t pb, int threads /* = MAX_THREADS */)
{
  std::vector<int32_t> primes;
  PrimeSieve ps;
  ps.generatePrimes(0, pb, &primes);
  int64_t sum = 0;

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for reduction(+: sum) schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = a + 1; i <= c; i++)
  {
    int64_t sum2 = 0;
    int64_t x2 = x / primes[i - 1];
    int64_t bi = pi_bsearch(primes.begin(), primes.end(), isqrt(x2));

    for (int64_t j = i; j <= bi; j++)
      sum2 -= pi_bsearch(primes.begin(), primes.end(), x2 / primes[j - 1]) - (j - 1);

    sum += sum2;
  }

  return sum;
}

} // namespace primecount
