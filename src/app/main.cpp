///
/// @file   main.cpp
/// @brief  primecount console application.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <ptypes.hpp>
#include <utils.hpp>

#include <stdexcept>
#include <iostream>
#include <limits>
#include <stdint.h>

using namespace std;
using namespace primecount;

namespace {

int64_t int64_cast(maxint_t x)
{
  if (x > numeric_limits<int64_t>::max())
    throw primecount_error("this is a 63-bit function, x must be < 2^63");
  return (int64_t) x;
}

} // namespace

int main (int argc, char* argv[])
{
  PrimeCountOptions pco = parseOptions(argc, argv);
  maxint_t res = 0;
  double time = get_wtime();

  try
  {
    switch (pco.option)
    {
      case OPTION_DELEGLISE_RIVAT:
        res = pi_deleglise_rivat(int64_cast(pco.x), pco.threads); break;
      case OPTION_DELEGLISE_RIVAT1:
        res = pi_deleglise_rivat1(int64_cast(pco.x)); break;
      case OPTION_DELEGLISE_RIVAT2:
        res = pi_deleglise_rivat2(int64_cast(pco.x)); break;
      case OPTION_DELEGLISE_RIVAT3:
        res = pi_deleglise_rivat3(int64_cast(pco.x)); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1:
        res = pi_deleglise_rivat_parallel1(int64_cast(pco.x), pco.threads); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL2:
        res = pi_deleglise_rivat_parallel2(int64_cast(pco.x), pco.threads); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL3:
        res = pi_deleglise_rivat_parallel3(int64_cast(pco.x), pco.threads); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(int64_cast(pco.x), pco.threads); break;
      case OPTION_LEHMER:
        res = pi_lehmer(int64_cast(pco.x), pco.threads); break;
      case OPTION_LEHMER2:
        res = pi_lehmer2(int64_cast(pco.x), pco.threads); break;
      case OPTION_LMO:
        res = pi_lmo(int64_cast(pco.x), pco.threads); break;
      case OPTION_LMO1:
        res = pi_lmo1(int64_cast(pco.x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(int64_cast(pco.x)); break;
      case OPTION_LMO3:
        res = pi_lmo3(int64_cast(pco.x)); break;
      case OPTION_LMO4:
        res = pi_lmo4(int64_cast(pco.x)); break;
      case OPTION_LMO5:
        res = pi_lmo5(int64_cast(pco.x)); break;
      case OPTION_LMO_PARALLEL1:
        res = pi_lmo_parallel1(int64_cast(pco.x), pco.threads); break;
      case OPTION_LMO_PARALLEL2:
        res = pi_lmo_parallel2(int64_cast(pco.x), pco.threads); break;
      case OPTION_LMO_PARALLEL3:
        res = pi_lmo_parallel3(int64_cast(pco.x), pco.threads); break;
      case OPTION_MEISSEL:
        res = pi_meissel(int64_cast(pco.x), pco.threads); break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(int64_cast(pco.x), pco.threads); break;
      case OPTION_PHI:
        res = phi(int64_cast(pco.x), int64_cast(pco.a)); break;
      case OPTION_PI:
        res = pi(int64_cast(pco.x), pco.threads); break;
      case OPTION_LI:
        res = Li(int64_cast(pco.x)); break;
      case OPTION_LIINV:
        res = Li_inverse(int64_cast(pco.x)); break;
      case OPTION_NTHPRIME:
        res = nth_prime(int64_cast(pco.x), pco.threads); break;
    }
  }
  catch (exception& e)
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  cout << res << endl;
  if (pco.time)
    cout << "Seconds: " << get_wtime() - time << endl;

  return 0;
}
