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
  {
    cout << str << endl;

    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      outfile << str << endl;
      outfile.close();
    }
  }
}

void print_log(const string& str)
{
  if (is_print())
  {
    cout << str << endl;
  }

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << str << endl;
    outfile.close();
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
}

void print_log(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
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

void print(const string& res_str, maxint_t res, double time)
{
  double seconds = get_wtime() - time;

  if (is_print())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << res_str << " = " << res << endl;
    print_seconds(seconds);

    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      outfile << "Status: 100%" << endl;
      outfile << res_str << " = " << res << endl;
      outfile << "Seconds: " << fixed << setprecision(3) << seconds << endl;
      outfile.close();
    }
  }
}

void print_seconds(double seconds)
{
  cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

} // namespace
