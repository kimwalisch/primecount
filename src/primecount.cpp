#include "cmdoptions.h"

#include <primecount.h>
#include <iostream>
#include <stdint.h>

using namespace primecount;

int main (int argc, char* argv[])
{
  // process command-line options
  PrimeCountSettings pcs = processOptions(argc, argv);
  int64_t res = 0;

  switch (pcs.option)
  {      
    case OPTION_PRIMESIEVE: res = pi_primesieve(pcs.x, pcs.threads); break;
    case OPTION_LEGENDRE:   res = pi_legendre  (pcs.x, pcs.threads); break;
    case OPTION_MEISSEL:    res = pi_meissel   (pcs.x, pcs.threads); break;
    case OPTION_LEHMER:     res = pi_lehmer    (pcs.x, pcs.threads); break;
    case OPTION_LI:         res = Li           (pcs.x); break;
    case OPTION_LIINV:      res = Li_inverse   (pcs.x); break;
    case OPTION_NTHPRIME:   res = nth_prime    (pcs.x, pcs.threads); break;
  }

  std::cout << res << std::endl;
  return 0;
}
