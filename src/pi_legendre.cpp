#include "utils/isqrt.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_legendre(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2) return 0;
  int64_t a = pi_primesieve(isqrt(x), threads);
  return a + phi(x, a, threads) - 1;
}

} // namespace primecount

#ifdef MAIN1

#include "utils/ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_legendre x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Legendre's prime counting formula." << std::endl
            << "The algorithm's complexity is O(x) operations and O(x^0.5) space."        << std::endl;
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

  std::cout << primecount::pi_legendre(x) << std::endl;
  return 0;
}

#endif /* MAIN */
