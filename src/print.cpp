///
/// @file  print.cpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
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

} // naespace

namespace primecount {

// backup.cpp
std::string backup_file();

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

void print_seconds(double seconds)
{
  cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

void print_status(double percent, maxint_t x)
{
  if (is_print())
    cout << "Status: " << fixed << setprecision(get_status_precision(x)) << percent << '%' << flush;
}

void print_resume(double percent, maxint_t x)
{
  print("Resuming from " + backup_file());

  // When the primecount.backup file is loaded from hard disk we
  // also verify the MD5 checksum and abort if there is a mismatch.
  // print_resume() is always called after that primecount.backup
  // has been loaded, hence we know the MD5 checksum is correct.
  print("MD5 checksum: OK");
  print_status(percent, x);
}

void print(const string& str)
{
  if (is_print())
    cout << str << endl;

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
    outfile << str << endl;
}

void print(const string& str, maxint_t res)
{
  if (is_print())
    cout << str << " = " << res << endl;

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
    outfile << str << " = " << res << endl;
}

void print(const string& str, maxint_t res, double time)
{
  if (is_print())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << str << " = " << res << endl;
    print_seconds(get_time() - time);
  }

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << "Status: 100%" << endl;
    outfile << str << " = " << res << endl;
    outfile << "Seconds: " << fixed << setprecision(3) << get_time() - time << endl;
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
  }

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
    }
  }
}

void print_gourdon(maxint_t x, int64_t y, int threads)
{
  if (print_variables())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }

  if (print_variables())
  {
    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
      outfile << "threads = " << threads << endl;
      outfile << endl;
    }
  }
}

void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads)
{
  if (print_variables())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "k = " << k << endl;
    cout << "x_star = " << get_x_star_gourdon(x, y) << endl;
    cout << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    cout << "alpha_z = " << fixed << setprecision(3) << get_alpha_z(y, z) << endl;
    print_threads(threads);
    cout << endl;
  }

  if (print_variables())
  {
    ofstream outfile("primecount.log", ofstream::out | ofstream::app);

    if (outfile.is_open())
    {
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "k = " << k << endl;
      outfile << "x_star = " << get_x_star_gourdon(x, y) << endl;
      outfile << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
      outfile << "alpha_z = " << fixed << setprecision(3) << get_alpha_z(y, z) << endl;
      outfile << "threads = " << threads << endl;
      outfile << endl;
    }
  }
}

void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, double alpha_y, double alpha_z, int threads)
{
  if (is_print())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "k = " << k << endl;
    cout << "x_star = " << get_x_star_gourdon(x, y) << endl;
    cout << "alpha_y = " << fixed << setprecision(3) << alpha_y << endl;
    cout << "alpha_z = " << fixed << setprecision(3) << alpha_z << endl;
    print_threads(threads);
  }

  ofstream outfile("primecount.log", ofstream::out | ofstream::app);

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "z = " << z << endl;
    outfile << "k = " << k << endl;
    outfile << "x_star = " << get_x_star_gourdon(x, y) << endl;
    outfile << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    outfile << "alpha_z = " << fixed << setprecision(3) << get_alpha_z(y, z) << endl;
    outfile << "threads = " << threads << endl;
  }
}

} // namespace
