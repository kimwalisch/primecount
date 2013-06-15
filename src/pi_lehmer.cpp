#include "isqrt.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_lehmer(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t c = pi_meissel(isqrt3(x));
  int64_t a = pi_meissel(isqrt4(x));
  int64_t b = pi_meissel(isqrt(x));

  int64_t sum = 0;

  sum += phi(x, a, threads);
  sum += P2(x, a, b, isqrt(x), threads);
  sum += P3(x, a, c, isqrt(x), threads);

  return sum;
}

} // namespace primecount
