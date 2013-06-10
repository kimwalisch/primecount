#include "../utils/ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_legendre x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Legendre's prime counting formula." << std::endl
            << "The complexity is O(x) operations and O(x^0.5) space."                    << std::endl;
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

  std::cout << legendre::pi(x) << std::endl;
  return 0;
}

