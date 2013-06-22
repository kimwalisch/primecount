#include "cmdoptions.h"

#include <primecount.h>
#include <iostream>
#include <stdint.h>

int main (int argc, char* argv[])
{
  // process command-line options
  PrimeCountSettings pcs = processOptions(argc, argv);
  int64_t res = 0;

  switch (pcs.method)
  {
    case 0: res = primecount::pi_primesieve(pcs.x, pcs.threads); break;
    case 1: res = primecount::pi_legendre(pcs.x, pcs.threads); break;
    case 2: res = primecount::pi_meissel(pcs.x, pcs.threads); break;
    case 3: res = primecount::pi_lehmer(pcs.x, pcs.threads); break;
    case 4: res = primecount::Li(pcs.x); break;
    case 5: res = primecount::Li_inverse(pcs.x); break;
    case 6: res = primecount::nth_prime(pcs.x, pcs.threads); break;
  }

  std::cout << res << std::endl;
  return 0;
}
