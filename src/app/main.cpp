///
/// @file   main.cpp
/// @brief  primecount console application.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"
#include "../deleglise-rivat/S2.hpp"

#include <primecount-internal.hpp>
#include <primecount.hpp>
#include <pmath.hpp>
#include <int128.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>

#include <stdint.h>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

using namespace std;
using namespace primecount;

namespace primecount {

int64_t int64_cast(maxint_t x)
{
  if (x > numeric_limits<int64_t>::max())
    throw primecount_error("this is a 63-bit function, x must be < 2^63");
  return (int64_t) x;
}

maxint_t P2(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("P2(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);

  if (x <= numeric_limits<int64_t>::max())
    return P2((int64_t) x, y, threads);
  else
    return P2(x, y, threads);
}

maxint_t S1(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("S1(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S1((int64_t) x, y, c, threads);
  else
    return S1(x, y, c, threads);
}

maxint_t S2_trivial(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("S2_trivial(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_trivial((int64_t) x, y, z, c, threads);
  else
    return S2_trivial(x, y, z, c, threads);
}

maxint_t S2_easy(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("S2_easy(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_easy((int64_t) x, y, z, c, threads);
  else
    return S2_easy(x, y, z, c, threads);
}

maxint_t S2_hard(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("S2_hard(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  // TODO: find better S2_hard approximation formula
  maxint_t s2_hard_approx = Li(x);

  if (x <= numeric_limits<int64_t>::max())
    return S2_hard((int64_t) x, y, z, c, (int64_t) s2_hard_approx, threads);
  else
    return S2_hard(x, y, z, c, s2_hard_approx, threads);
}

} // namespace primecount

int main (int argc, char* argv[])
{
  PrimeCountOptions pco = parseOptions(argc, argv);
  double time = get_wtime();

  maxint_t x = pco.x;
  maxint_t res = 0;
  int threads = pco.threads;

  try
  {
    switch (pco.option)
    {
      case OPTION_DELEGLISE_RIVAT:
        res = pi_deleglise_rivat(x, threads); break;
      case OPTION_DELEGLISE_RIVAT1:
        res = pi_deleglise_rivat1(int64_cast(x)); break;
      case OPTION_DELEGLISE_RIVAT2:
        res = pi_deleglise_rivat2(int64_cast(x)); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1:
        res = pi_deleglise_rivat_parallel1(int64_cast(x), threads); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL2:
        res = pi_deleglise_rivat_parallel2(int64_cast(x), threads); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(int64_cast(x), threads); break;
      case OPTION_LEHMER:
        res = pi_lehmer(int64_cast(x), threads); break;
      case OPTION_LMO:
        res = pi_lmo(int64_cast(x), threads); break;
      case OPTION_LMO1:
        res = pi_lmo1(int64_cast(x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(int64_cast(x)); break;
      case OPTION_LMO3:
        res = pi_lmo3(int64_cast(x)); break;
      case OPTION_LMO4:
        res = pi_lmo4(int64_cast(x)); break;
      case OPTION_LMO5:
        res = pi_lmo5(int64_cast(x)); break;
      case OPTION_LMO_PARALLEL1:
        res = pi_lmo_parallel1(int64_cast(x), threads); break;
      case OPTION_LMO_PARALLEL2:
        res = pi_lmo_parallel2(int64_cast(x), threads); break;
      case OPTION_LMO_PARALLEL3:
        res = pi_lmo_parallel3(int64_cast(x), threads); break;
      case OPTION_MEISSEL:
        res = pi_meissel(int64_cast(x), threads); break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(int64_cast(x), threads); break;
      case OPTION_P2:
        res = P2(x, threads); break;
      case OPTION_PI:
        res = pi(x, threads); break;
      case OPTION_LI:
        res = Li(int64_cast(x)); break;
      case OPTION_LIINV:
        res = Li_inverse(int64_cast(x)); break;
      case OPTION_NTHPRIME:
        res = nth_prime(int64_cast(x), threads); break;
      case OPTION_S1:
        res = S1(x, threads); break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads); break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads); break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads); break;

#ifdef HAVE_INT128_T
      case OPTION_DELEGLISE_RIVAT_PARALLEL3:
        res = pi_deleglise_rivat_parallel3(x, threads); break;
#endif
    }
  }
  catch (bad_alloc&)
  {
    cerr << "Error: failed to allocate memory, your system most likely does" << endl
         << "       not have enough memory to run this computation." << endl;
    return 1;
  }
  catch (exception& e)
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  if (print_result())
  {
    if (print_status())
      cout << endl;
    cout << res << endl;
    if (pco.time)
      print_seconds(get_wtime() - time);
  }

  return 0;
}
