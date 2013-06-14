#include "utils/isqrt.h"
#include "utils/PrimeCountVector.h"

#include <primesieve/soe/PrimeSieve.h>
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
  int64_t pix = 0;
  int64_t old = 0;
  int64_t sum = phi(x, a, threads) + ((b + a - 2) * (b - a + 1) / 2);

  PrimeCountVector<uint32_t> primes;
  PrimeSieve ps;
  ps.generatePrimes(0, isqrt(x), &primes);
  
  // for use as array indexes
  a--; b--; c--;

  for (int64_t i = b; i > a; i--)
  {
    int64_t x2 = x / primes[i];
    pix += ps.countPrimes(old + 1, x2);
    old = x2;
    sum -= pix;

    if (i <= c) {
      int64_t bi = primes.pi( isqrt(x2) ) - 1;
      for (int64_t j = i; j <= bi; j++)
        sum = sum - primes.pi( x2 / primes[j]) + (j + 1) - 1;
    }
  }

  return sum;
}

} // namespace primecount

#ifdef MAIN3

#include "utils/ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_lehmer x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Lehmer0's prime counting formula." << std::endl
            << "The complexity is O( TODO ) operations and O( TODO ) space."     << std::endl;
  std::exit(1);
}

} // end namespace

int main (int argc, char* argv[])
{
  if (argc != 2)
    help();

  ExpressionParser<int64_t> parser;
  int64_t x = 0;
  try {
    x = parser.eval(argv[1]);
  }
  catch (parser_error&) {
    help();
  }

  std::cout << primecount::pi_lehmer(x) << std::endl;
  return 0;
}

#endif /* MAIN */
