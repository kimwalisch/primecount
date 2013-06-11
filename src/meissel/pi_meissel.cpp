#include "../utils/isqrt.h"

#include <primecount.h>

namespace primecount {

int64_t pi_meissel(int64_t x, int threads /* = MAX_THREADS */)
{
  int64_t sum = 0;
  if (x >= 2)
  {
    int64_t a_sqrt3 = isqrt3(x);
    int64_t b_sqrt2 = isqrt(x);
    int64_t a = pi_legendre(a_sqrt3);

    sum = phi(x, a, threads) + P2(x, a, a_sqrt3, b_sqrt2, threads);
  }
  return sum;
}

} // namespace primecount
