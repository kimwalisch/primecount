#include "isqrt.h"

#include <primecount.h>
#include <primesieve/soe/ParallelPrimeSieve.h>
#include <stdint.h>

namespace primecount {

/// This function calculates the nth prime using a combination of an
/// efficient implementation of the prime counting function pi(x)
/// (currently Lehmer's method) and an implementation of the segmented
/// sieve of Eratosthenes (the author's primesieve library).
///
int64_t nth_prime(int64_t n, int threads /* = MAX_THREADS */)
{
  if (n < 1)
    return 0;

  int64_t prime_approx = 0;
  int64_t count_approx = 0;

  if (n > 10000)
  {
    // this is < prime[n] up to ~ 10^316
    prime_approx = Li_inverse(n);
    count_approx = pi(prime_approx, threads);
  }

  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS)
    pps.setNumThreads( threads );
  int64_t prime = pps.nthPrime(prime_approx, n - count_approx);

  return prime;
}

} // namespace primecount
