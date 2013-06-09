#include "../meissel/pi.h"
#include "../utils/ExpressionParser.h"

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

  std::cout << meissel::pi(x) << std::endl;
  return 0;
}

