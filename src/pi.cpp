#include <primecount.h>
#include <stdint.h>

namespace primecount {

/// Alias for fastest pi(x) implementation
int64_t pi(int64_t x, int threads /* = MAX_THREADS */)
{
  return pi_lehmer(x, threads);
}

} // namespace primecount
