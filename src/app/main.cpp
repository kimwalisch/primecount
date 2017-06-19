///
/// @file   main.cpp
/// @brief  primecount console application
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2.hpp>
#include <json.hpp>

#include <stdint.h>
#include <exception>
#include <fstream>
#include <iomanip>
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

maxint_t P2(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primecount_error("P2(x): x must be <= " + limit);

  if (is_print())
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

  if (is_print())
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

  if (is_print())
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

  if (is_print())
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

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_hard((int64_t) x, y, z, c, (int64_t) Ri(x), threads);
  else
    return S2_hard(x, y, z, c, Ri(x), threads);
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

    if (!opt.resume)
      backup_command(argc, argv);
    else
    {
      vector<string> args = get_backup_command();
      vector<char*> cargs;

      for (string& s : args)
        cargs.push_back((char*) s.c_str());

      opt = parseOptions((int) cargs.size(), cargs.data());
    }

    double time = get_wtime();

    auto x = opt.x;
    auto a = opt.a;
    auto threads = opt.threads;
    maxint_t res = 0;

    switch (opt.option)
    {
      case OPTION_DELEGLISE_RIVAT:
      {
#ifdef HAVE_INT128_T
        if (x > numeric_limits<int64_t>::max())
          res = pi_deleglise_rivat_parallel2(x, threads, &time);
        else
#endif
        res = pi_deleglise_rivat_parallel1(to_int64(x), threads, &time);
        break;
      }
      case OPTION_DELEGLISE_RIVAT1:
        res = pi_deleglise_rivat1(to_int64(x)); break;
      case OPTION_DELEGLISE_RIVAT2:
        res = pi_deleglise_rivat2(to_int64(x)); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1:
        res = pi_deleglise_rivat_parallel1(to_int64(x), threads, &time); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(to_int64(x), threads); break;
      case OPTION_LEHMER:
        res = pi_lehmer(to_int64(x), threads); break;
      case OPTION_LMO:
        res = pi_lmo(to_int64(x), threads); break;
      case OPTION_LMO1:
        res = pi_lmo1(to_int64(x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(to_int64(x)); break;
      case OPTION_LMO3:
        res = pi_lmo3(to_int64(x)); break;
      case OPTION_LMO4:
        res = pi_lmo4(to_int64(x)); break;
      case OPTION_LMO5:
        res = pi_lmo5(to_int64(x)); break;
      case OPTION_LMO_PARALLEL:
        res = pi_lmo_parallel(to_int64(x), threads); break;
      case OPTION_MEISSEL:
        res = pi_meissel(to_int64(x), threads); break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(to_int64(x), threads); break;
      case OPTION_P2:
        res = P2(x, threads); break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); break;
      case OPTION_PI:
        res = pi(x, threads); break;
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
      case OPTION_S1:
        res = S1(x, threads); break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads); break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads); break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads); break;
#ifdef HAVE_INT128_T
      case OPTION_DELEGLISE_RIVAT_PARALLEL2:
        res = pi_deleglise_rivat_parallel2(x, threads, &time); break;
#endif
    }

    if (print_result())
    {
      ofstream outfile("primecount.log", ofstream::out | ofstream::app);
      double seconds = get_wtime() - time;

      if (is_print())
        cout << endl;

      cout << res << endl;
      
      if (opt.time)
        cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
  
      outfile << endl << res << endl;
      outfile << "Seconds: " << fixed << setprecision(3) << seconds << endl;
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
