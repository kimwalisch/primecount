#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <vector>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

int64_t P2(int64_t x, int64_t a, int64_t b, int64_t pb, int threads /* = MAX_THREADS */)
{
  std::vector<int32_t> primes;
  std::vector<int64_t> counts;
  PrimeSieve ps;
  ps.generatePrimes(0, pb, &primes);
  counts.resize( primes.size() );

  int64_t size = static_cast<int64_t>( primes.size() );

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for private(ps) schedule(dynamic) num_threads(threads)
#endif
  for (int64_t i = size; i > a; i--)
  {
    int64_t start = (i < size) ? x / primes[i] : 0;
    int64_t stop = x / primes[i - 1];
    counts[i - 1] = ps.countPrimes(start + 1, stop);
  }

  int64_t sum = (b + a - 2) * (b - a + 1) / 2;
  int64_t pix = 0;

  for (int64_t i = size; i > a; i--)
  {
    pix += counts[i - 1];
    sum -= pix;
  }

  return sum;
}

} // namespace primecount
