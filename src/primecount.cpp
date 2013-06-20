#include "cmdoptions.h"

#include <primecount.h>
#include <iostream>
#include <stdint.h>

int main (int argc, char* argv[])
{
  // process command-line options
  PrimeCountSettings pcs = processOptions(argc, argv);
  int64_t pix = 0;

  switch (pcs.method)
  {
    case 0: pix = primecount::pi_primesieve(pcs.x, pcs.threads); break;
    case 1: pix = primecount::pi_legendre(pcs.x, pcs.threads); break;
    case 2: pix = primecount::pi_meissel(pcs.x, pcs.threads); break;
    case 3: pix = primecount::pi_lehmer(pcs.x, pcs.threads); break;
  }

  std::cout << pix << std::endl;
  return 0;
}
