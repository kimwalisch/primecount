#include "isqrt.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_legendre(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t a = pi_primesieve(isqrt(x), threads);
  int64_t sum = a + phi(x, a, threads) - 1;

  return sum;
}

} // namespace primecount
