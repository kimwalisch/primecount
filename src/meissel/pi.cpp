#include "../utils/isqrt.h"

#include <primecount.h>

namespace meissel {

int64_t pi(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2) return 0;
  int64_t a = legendre::pi(isqrt3(x));
  return legendre::phi(x, a, threads) + Pkxa::P2(x, a, threads);
}

} // namespace meissel
