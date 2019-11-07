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

string basename(string path)
{
  while (path.back() == '/' ||
         path.back() == '\\')
  {
    path.pop_back();
  }

  size_t pos = path.find_last_of("/\\");

  if (pos != string::npos)
    return path.substr(pos + 1);

  return path;
}

void result_txt(int argc,
                char* argv[],
                maxint_t res,
                int threads,
                double seconds)
{
  ofstream resfile("results.txt", ofstream::out | ofstream::app);

  if (resfile.is_open())
  {
    auto json = load_backup();

    // Don't put primecount --resume into results file
    if (json.count("command") > 0)
    {
      int i = 0;
      int size = 0;

      for (auto& c : json["command"])
        size += 1;

      for (auto& c : json["command"])
      {
        resfile << c.get<std::string>();
        i += 1;
        if (i < size)
          resfile << " ";
      }
    }
    else
    {
      resfile << basename(argv[0]);
      for (int i = 1; i < argc; i++)
        resfile << " " << argv[i];
    }

    resfile << endl;
    resfile << "Result: " << res << endl;
    resfile << "Threads: " << threads << endl;
    resfile << "Seconds: " << fixed << setprecision(3) << seconds << endl;
    resfile << "Date: " << dateTime() << endl << endl;
  }
}

void log_footer()
{
  ofstream logfile("primecount.log", ofstream::out | ofstream::app);

  if (logfile.is_open())
    logfile << endl << string(60, '=') << endl;
}

void log_result(maxint_t res, double seconds)
{
  ofstream logfile("primecount.log", ofstream::out | ofstream::app);

  if (logfile.is_open())
  {
    logfile << endl << res << endl;
    logfile << "Seconds: " << fixed << setprecision(3) << seconds << endl;
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

    auto x = opt.x;
    auto a = opt.a;
    auto time = get_time();
    auto threads = get_num_threads();
    maxint_t res = 0;

    switch (opt.option)
    {
      case OPTION_DEFAULT:
        res = pi(x, threads);
        time = backup_time(time); break;
      case OPTION_GOURDON:
        res = pi_gourdon(x, threads);
        time = backup_time(time); break;
      case OPTION_GOURDON_64:
        res = pi_gourdon_64(to_int64(x), threads);
        time = backup_time(time); break;
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
        res = nth_prime(to_int64(x), threads);
        time = backup_time(time); break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); break;
      case OPTION_AC:
        res = AC(x, threads);
        time = backup_time(time, "AC"); break;
      case OPTION_B:
        res = B(x, threads);
        time = backup_time(time, "B"); break;
      case OPTION_D:
        res = D(x, threads);
        time = backup_time(time, "D"); break;
      case OPTION_PHI0:
        res = Phi0(x, threads);
        time = backup_time(time, "Phi0"); break;
      case OPTION_SIGMA:
        res = Sigma(x, threads);
        time = backup_time(time, "Sigma"); break;
#ifdef HAVE_INT128_T
      case OPTION_GOURDON_128:
        res = pi_gourdon_128(x, threads);
        time = backup_time(time); break;
#endif
    }

    double seconds = get_time() - time;
    result_txt(argc, argv, res, threads, seconds);

    if (print_result())
    {
      log_result(res, seconds);

      if (is_print())
        cout << endl;

      cout << res << endl;

      if (opt.time)
        cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
    }

    log_footer();
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
