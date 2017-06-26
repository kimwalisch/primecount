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
#include <stdio.h>
#include <time.h>

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

namespace {

// Get current date & time
// http://www.cplusplus.com/reference/ctime/strftime
const string dateTime()
{
  time_t rawtime;
  struct tm* timeinfo;
  char buffer[80];

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, sizeof(buffer), "%Y-%m-%d %X", timeinfo);

  return buffer;
}

void result_log(maxint_t res, double seconds)
{
  ofstream logfile("primecount.log", ofstream::out | ofstream::app);

  if (logfile.is_open())
  {
    logfile << endl << res << endl;
    logfile << "Seconds: " << fixed << setprecision(3) << seconds << endl;
  }
}

void result_txt(string& formula,
                maxint_t res,
                CmdOptions& opt,
                double seconds)
{
  ofstream resfile("results.txt", ofstream::out | ofstream::app);

  if (resfile.is_open())
  {
    resfile << formula << " = " << res << endl;
    resfile << "Command: primecount " << opt.args << endl;
    resfile << "Threads: " << opt.threads << endl;
    resfile << "Seconds: " << seconds << endl;
    resfile << "Date: " << dateTime() << endl << endl;
  }
}

} // namespace

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

    if (!opt.is_resume())
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
    string formula;

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
        formula = "pi(" + opt.x_str + ")";
        break;
      }
      case OPTION_DELEGLISE_RIVAT1:
        res = pi_deleglise_rivat1(to_int64(x)); formula = "pi(" + opt.x_str + ")"; break;
      case OPTION_DELEGLISE_RIVAT2:
        res = pi_deleglise_rivat2(to_int64(x)); formula = "pi(" + opt.x_str + ")"; break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1:
        res = pi_deleglise_rivat_parallel1(to_int64(x), threads, &time); formula = "pi(" + opt.x_str + ")"; break;
      case OPTION_LEGENDRE:
        res = pi_legendre(to_int64(x), threads); formula = "pi_legendre(" + opt.x_str + ")"; break;
      case OPTION_LEHMER:
        res = pi_lehmer(to_int64(x), threads); formula = "pi_lehmer(" + opt.x_str + ")"; break;
      case OPTION_LMO:
        res = pi_lmo(to_int64(x), threads); formula = "pi_lmo(" + opt.x_str + ")"; break;
      case OPTION_LMO1:
        res = pi_lmo1(to_int64(x)); formula = "pi_lmo1(" + opt.x_str + ")"; break;
      case OPTION_LMO2:
        res = pi_lmo2(to_int64(x)); formula = "pi_lmo2(" + opt.x_str + ")"; break;
      case OPTION_LMO3:
        res = pi_lmo3(to_int64(x)); formula = "pi_lmo3(" + opt.x_str + ")"; break;
      case OPTION_LMO4:
        res = pi_lmo4(to_int64(x)); formula = "pi_lmo4(" + opt.x_str + ")"; break;
      case OPTION_LMO5:
        res = pi_lmo5(to_int64(x)); formula = "pi_lmo5(" + opt.x_str + ")"; break;
      case OPTION_LMO_PARALLEL:
        res = pi_lmo_parallel(to_int64(x), threads); formula = "pi_lmo(" + opt.x_str + ")"; break;
      case OPTION_MEISSEL:
        res = pi_meissel(to_int64(x), threads); formula = "pi_meissel(" + opt.x_str + ")"; break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(to_int64(x), threads); formula = "pi_primesieve(" + opt.x_str + ")"; break;
      case OPTION_P2:
        res = P2(x, threads); formula = "P2(" + opt.x_str + ")"; break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); formula = "phi(" + opt.x_str + ", " + to_string(a) + ")"; break;
      case OPTION_PI:
        res = pi(x, threads); formula = "pi(" + opt.x_str + ")"; break;
      case OPTION_LI:
        res = Li(x); formula = "Li(" + opt.x_str + ")"; break;
      case OPTION_LIINV:
        res = Li_inverse(x); formula = "Li_inverse(" + opt.x_str + ")"; break;
      case OPTION_RI:
        res = Ri(x); formula = "Ri(" + opt.x_str + ")"; break;
      case OPTION_RIINV:
        res = Ri_inverse(x); formula = "Ri_inverse(" + opt.x_str + ")"; break;
      case OPTION_NTHPRIME:
        res = nth_prime(to_int64(x), threads); formula = "nth_prime(" + opt.x_str + ")"; break;
      case OPTION_S1:
        res = S1(x, threads); formula = "S1(" + opt.x_str + ")"; break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads); formula = "S2_easy(" + opt.x_str + ")"; break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads); formula = "S2_hard(" + opt.x_str + ")"; break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads); formula = "S2_trivial(" + opt.x_str + ")"; break;
#ifdef HAVE_INT128_T
      case OPTION_DELEGLISE_RIVAT_PARALLEL2:
        res = pi_deleglise_rivat_parallel2(x, threads, &time); formula = "pi(" + opt.x_str + ")"; break;
#endif
    }

    double seconds = get_wtime() - time;
    result_txt(formula, res, opt, seconds);

    if (print_result())
    {
      if (is_print())
        cout << endl;

      cout << res << endl;

      if (opt.time)
        cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;

      result_log(res, seconds);
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
