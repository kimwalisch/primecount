#include "../utils/isqrt.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>

namespace primecount {

int64_t pi_legendre(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2) return 0;
  PrimeSieve ps;
  int64_t a = ps.countPrimes(0, isqrt(x));
  return a + phi(x, a, threads) - 1;
}

} // namespace primecount
