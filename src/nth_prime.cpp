#include "isqrt.h"

#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>
#include <primesieve/soe/PrimeSieveCallback.h>
#include <primesieve/soe/stop_primesieve.h>
#include <stdint.h>

namespace {

/// This class is used to generate n primes and
/// then stop by throwing an exception.
///
class FindNthPrime : public PrimeSieveCallback<uint64_t> {
public:
  FindNthPrime(uint64_t n) : n_(n), nthPrime_(0)
  { }
  void callback(uint64_t prime)
  {
    if (--n_ == 0)
    {
      nthPrime_ = prime;
      throw stop_primesieve();
    }
  }
  int64_t getNthPrime() const
  {
    return nthPrime_;
  }
private:
  uint64_t n_;
  uint64_t nthPrime_;
};

} // namespace

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

  int64_t nth_prime_approximation = Li_inverse(n);
  int64_t a = nth_prime_approximation;
  int64_t pia = pi(a, threads);

  FindNthPrime findNthPrime(n - pia);
  try {
    int64_t start = a;
    int64_t stop = a + isqrt(a) * 2 + 10000;
    PrimeSieve ps;
    ps.generatePrimes(start, stop, &findNthPrime);
  }
  catch (stop_primesieve&) {
    return findNthPrime.getNthPrime();
  }

  return -1;
}

} // namespace primecount
