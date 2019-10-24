///
/// @file   main.cpp
/// @brief  primecount console application
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <PhiTiny.hpp>
#include <print.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

#ifdef HAVE_MPI
  #include <mpi.h>
#endif

using namespace std;
using namespace primecount;

namespace primecount {

int64_t to_int64(maxint_t x)
{
  if (x > numeric_limits<int64_t>::max())
    throw primecount_error("x must be < 2^63");
  return (int64_t) x;
}

maxint_t AC(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("AC(x): x must be <= " + to_str(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= numeric_limits<int64_t>::max())
    return AC((int64_t) x, y, z, k, threads);
  else
    return AC(x, y, z, k, threads);
}

maxint_t B(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("B(x): x must be <= " + to_str(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= numeric_limits<int64_t>::max())
    return B((int64_t) x, y, threads);
  else
    return B(x, y, threads);
}

maxint_t D(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("D(x): x must be <= " + to_str(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= numeric_limits<int64_t>::max())
    return D((int64_t) x, y, z, k, (int64_t) Ri(x), threads);
  else
    return D(x, y, z, k, Ri(x), threads);
}

maxint_t Phi0(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("Phi0(x): x must be <= " + to_str(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= numeric_limits<int64_t>::max())
    return Phi0((int64_t) x, y, z, k, threads);
  else
    return Phi0(x, y, z, k, threads);
}

maxint_t Sigma(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("Sigma(x): x must be <= " + to_str(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= numeric_limits<int64_t>::max())
    return Sigma((int64_t) x, y, threads);
  else
    return Sigma(x, y, threads);
}

} // namespace

int main (int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  try
  {
    CmdOptions opt = parseOptions(argc, argv);
    double time = get_time();

    auto x = opt.x;
    auto a = opt.a;
    auto threads = get_num_threads();
    maxint_t res = 0;

    switch (opt.option)
    {
      case OPTION_DEFAULT:
        res = pi(x, threads); break;
      case OPTION_GOURDON:
        res = pi_gourdon(x, threads); break;
      case OPTION_GOURDON_64:
        res = pi_gourdon_64(to_int64(x), threads); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(to_int64(x), threads); break;
      case OPTION_MEISSEL:
        res = pi_meissel(to_int64(x), threads); break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(to_int64(x)); break;
      case OPTION_LI:
        res = Li(x); break;
      case OPTION_LIINV:
        res = Li_inverse(x); break;
      case OPTION_RI:
        res = Ri(x); break;
      case OPTION_RIINV:
        res = Ri_inverse(x); break;
      case OPTION_NTHPRIME:
        res = nth_prime(to_int64(x), threads); break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); break;
      case OPTION_AC:
        res = AC(x, threads); break;
      case OPTION_B:
        res = B(x, threads); break;
      case OPTION_D:
        res = D(x, threads); break;
      case OPTION_PHI0:
        res = Phi0(x, threads); break;
      case OPTION_SIGMA:
        res = Sigma(x, threads); break;
#ifdef HAVE_INT128_T
      case OPTION_GOURDON_128:
        res = pi_gourdon_128(x, threads); break;
#endif
    }

    if (print_result())
    {
      if (is_print())
        cout << endl;

      cout << res << endl;

      if (opt.time)
        print_seconds(get_time() - time);
    }
  }
  catch (exception& e)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    cerr << "primecount: " << e.what() << endl
         << "Try 'primecount --help' for more information." << endl;
    return 1;
  }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

  return 0;
}
