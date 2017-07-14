///
/// @file  print.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <backup.hpp>
#include <stdint.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_ = false;

bool print_variables_ = false;

}

namespace primecount {

void set_print(bool print)
{
#ifdef HAVE_MPI
  print_ = print && is_mpi_master_proc();
#else
  print_ = print;
#endif
}

void set_print_variables(bool print_variables)
{
#ifdef HAVE_MPI
  print_variables_ = print_variables && is_mpi_master_proc();
#else
  print_variables_ = print_variables;
#endif
}

bool print_result()
{
#ifdef HAVE_MPI
  return !print_variables() && is_mpi_master_proc();
#else
  return !print_variables();
#endif
}

void print_threads(int threads)
{
#ifdef HAVE_MPI
  cout << "processes = " << mpi_num_procs() << endl;
  cout << "threads = " << mpi_num_procs() << " * " << threads << endl;
#else
  cout << "threads = " << threads << endl;
#endif
}

bool is_print()
{
  return print_;
}

bool print_variables()
{
  return print_variables_;
}

void print(const string& str)
{
  if (is_print())
    cout << str << endl;
}

void print(maxint_t x, int64_t y, int threads)
{
  if (print_variables())
  {
    maxint_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }
}

void print(maxint_t x, int64_t y, int64_t c, int threads)
{
  if (print_variables())
  {
    maxint_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }
}

void print(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
{
  if (is_print())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    print_threads(threads);
  }
}

void print(const string& res_str, maxint_t res)
{
  if (is_print())
    cout << res_str << " = " << res << endl;
}

void print(const string& res_str, maxint_t res, double time)
{
  if (is_print())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << res_str << " = " << res << endl;
    cout << "Seconds: " << fixed << setprecision(3) << get_wtime() - time << endl;
  }
}

void print_log(const string& str)
{
  print(str);

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << str << endl;
    outfile.close();
  }
}

void print_log(maxint_t x, int64_t y, int threads)
{
  print(x, y, threads);

  if (print_variables())
  {
    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      maxint_t z = x / y;
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
      outfile << "threads = " << threads << endl;
      outfile << endl;
      outfile.close();
    }
  }
}

void print_log(maxint_t x, int64_t y, int64_t c, int threads)
{
  print(x, y, c, threads);

  if (print_variables())
  {
    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      maxint_t z = x / y;
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "c = " << c << endl;
      outfile << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
      outfile << "threads = " << threads << endl;
      outfile << endl;
      outfile.close();
    }
  }
}

void print_log(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
{
  print(x, y, z, c, alpha, threads);

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "z = " << z << endl;
    outfile << "c = " << c << endl;
    outfile << "alpha = " << fixed << setprecision(3) << alpha << endl;
    outfile << "threads = " << threads << endl;
    outfile.close();
  }
}

void print_log(const string& res_str, maxint_t res)
{
  print(res_str, res);

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << res_str << " = " << res << endl;
    outfile.close();
  }
}

void print_log(const string& res_str, maxint_t res, double time)
{
  print(res_str, res, time);

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << "Status: 100%" << endl;
    outfile << res_str << " = " << res << endl;
    outfile << "Seconds: " << fixed << setprecision(3) << get_wtime() - time << endl;
    outfile.close();
  }
}

void print_seconds(double seconds)
{
  if (is_print())
    cout << "Seconds: " << fixed << setprecision(3) << seconds << endl << endl;
}

void print_log_seconds(double seconds)
{
  print_seconds(seconds);

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << "Seconds: " << fixed << setprecision(3) << seconds << endl << endl;
    outfile.close();
  }
}

void print_status(double percent, maxint_t x)
{
  if (is_print())
    cout << "Status: " << fixed << setprecision(get_status_precision(x)) << percent << '%' << flush;
}

void print_resume(double percent, maxint_t x)
{
  print_log("Resuming from " + backup_file());
  print_status(percent, x);
}

} // namespace
