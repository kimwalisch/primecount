#include "utils/isqrt.h"

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t pi_meissel(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  int64_t pa_max = isqrt3(x);
  int64_t pb_max = isqrt(x);
  int64_t a = pi_legendre(pa_max);
  int64_t b = pi_legendre(pb_max);

  return phi(x, a, threads) + P2(x, a, b, pb_max, threads);
}

} // namespace primecount

#ifdef MAIN

#include "utils/ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_meissel x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Meissel's prime counting formula." << std::endl
            << "The complexity is O(x/(ln x)^3) operations and O(x^0.5/ln x) space."     << std::endl;
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

  std::cout << primecount::pi_meissel(x) << std::endl;
  return 0;
}

#endif /* MAIN */
