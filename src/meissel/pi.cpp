#include "pi.h"
#include "../utils/isqrt.h"
#include "../utils/PrimeSieveVector.h"
#include "../legendre/pi.h"
#include "../legendre/phi.h"

#include <primesieve/soe/PrimeSieve.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace meissel {

int64_t P2xa(int64_t b, int64_t c, int64_t x, const std::vector<uint32_t>& primes, int threads)
{
  int64_t pix = 0;
  int64_t old = 0;
  int64_t sum = 0;
  PrimeSieve ps;
  int size = static_cast<int>(primes.size());

#ifdef _OPENMP
  if (threads == MAX_THREADS)
    threads = omp_get_max_threads();
  #pragma omp parallel for private(ps) firstprivate(pix, old) reduction(+: sum) \
      num_threads(threads) schedule(dynamic)
#endif
  for (int i = size - 1; i >= 0; i--)
  {
    int64_t x2 = x / primes[i];
    if (old < x2)
      pix += ps.countPrimes(old + 1, x2);
    old = x2;
    sum += pix;
  }
  return (b + c - 2) * (b - c + 1) / 2 - sum;
}

int64_t pi(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t sqrt2 = isqrt(x);
  int64_t sqrt3 = isqrt3(x);

  PrimeSieveVector<uint32_t> primes;
  PrimeSieve ps;
  if (sqrt3 < sqrt2)
    ps.generatePrimes(sqrt3 + 1, sqrt2, &primes);

  int64_t b = legendre::pi(sqrt2);
  int64_t c = legendre::pi(sqrt3);

  return legendre::phi(x, c, threads) + P2xa(b, c, x, primes, threads);
}

} // namespace meissel
