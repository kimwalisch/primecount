///
/// @file  primecount.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <exception>
#include <iostream>
#include <stdint.h>

using namespace std;
using namespace primecount;

int main (int argc, char* argv[])
{
  PrimeCountOptions pco = parseOptions(argc, argv);
  int64_t res = 0;

  try
  {
    switch (pco.option)
    {
      case OPTION_DELEGLISE_RIVAT:           res = pi_deleglise_rivat          (pco.x, pco.threads); break;
      case OPTION_DELEGLISE_RIVAT1:          res = pi_deleglise_rivat1         (pco.x); break;
      case OPTION_DELEGLISE_RIVAT2:          res = pi_deleglise_rivat2         (pco.x); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1: res = pi_deleglise_rivat_parallel1(pco.x, pco.threads); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL2: res = pi_deleglise_rivat_parallel2(pco.x, pco.threads); break;
      case OPTION_LEGENDRE:                  res = pi_legendre                 (pco.x, pco.threads); break;
      case OPTION_LEHMER:                    res = pi_lehmer                   (pco.x, pco.threads); break;
      case OPTION_LEHMER2:                   res = pi_lehmer2                  (pco.x, pco.threads); break;
      case OPTION_LMO:                       res = pi_lmo                      (pco.x, pco.threads); break;
      case OPTION_LMO1:                      res = pi_lmo1                     (pco.x); break;
      case OPTION_LMO2:                      res = pi_lmo2                     (pco.x); break;
      case OPTION_LMO3:                      res = pi_lmo3                     (pco.x); break;
      case OPTION_LMO4:                      res = pi_lmo4                     (pco.x); break;
      case OPTION_LMO5:                      res = pi_lmo5                     (pco.x); break;
      case OPTION_LMO6:                      res = pi_lmo6                     (pco.x); break;
      case OPTION_LMO_PARALLEL1:             res = pi_lmo_parallel1            (pco.x, pco.threads); break;
      case OPTION_LMO_PARALLEL2:             res = pi_lmo_parallel2            (pco.x, pco.threads); break;
      case OPTION_LMO_PARALLEL3:             res = pi_lmo_parallel3            (pco.x, pco.threads); break;
      case OPTION_LMO_PARALLEL4:             res = pi_lmo_parallel4            (pco.x, pco.threads); break;
      case OPTION_MEISSEL:                   res = pi_meissel                  (pco.x, pco.threads); break;
      case OPTION_PRIMESIEVE:                res = pi_primesieve               (pco.x, pco.threads); break;
      case OPTION_PHI:                       res = phi                         (pco.x, pco.a); break;
      case OPTION_PI:                        res = pi                          (pco.x, pco.threads); break;
      case OPTION_LI:                        res = Li                          (pco.x); break;
      case OPTION_LIINV:                     res = Li_inverse                  (pco.x); break;
      case OPTION_NTHPRIME:                  res = nth_prime                   (pco.x, pco.threads); break;
    }
  }
  catch (exception& e)
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  cout << res << endl;
  return 0;
}
