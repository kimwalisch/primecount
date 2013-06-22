#include "P2.h"
#include "isqrt.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_meissel(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t a = pi_legendre(isqrt3(x));
  int64_t b = pi_legendre(isqrt(x));

  int64_t sum = 0;

  sum += phi(x, a, threads);
  sum += P2(x, a, b, isqrt(x), threads);

  return sum;
}

} // namespace primecount
