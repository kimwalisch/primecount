#include "ExpressionParser.h"

#include <primecount.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

namespace {

void help()
{
  std::cerr << "Usage: pi_lehmer x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Lehmer's prime counting formula." << std::endl
            << "The computational complexity is O(x/ln(x)^3) operations and"            << std::endl
            << "O(x^(1/3)/ln(x)) space."                                                << std::endl;
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
    primecount::test();
  }

  std::cout << primecount::pi(x) << std::endl;
  return 0;
}
