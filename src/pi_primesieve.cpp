#include "utils/isqrt.h"

#include <primecount.h>
#include <primesieve/soe/ParallelPrimeSieve.h>
#include <stdint.h>

namespace primecount {

int64_t pi_primesieve(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2) return 0;
  ParallelPrimeSieve pps;
  if (threads != MAX_THREADS) pps.setNumThreads(threads);
  return pps.countPrimes(0, x);
}

} // namespace primecount

#ifdef MAIN2

#include "utils/ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_primesieve x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using primesieve (multi-threaded segmented" << std::endl
            << "sieve of Eratothenes implementation). "                                     << std::endl
            << "The algorithm's complexity is O(ln ln x) operations and O(x^0.5) space."    << std::endl;
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

  std::cout << primecount::pi_primesieve(x) << std::endl;
  return 0;
}

#endif /* MAIN */
