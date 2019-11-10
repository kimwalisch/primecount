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

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_ = false;
bool print_variables_ = false;

bool is_print_variables()
{
  return print_variables_;
}

} // naespace

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

bool is_print()
{
  return print_;
}

/// The final combined result is always shown at
/// the end even if is_print = false. It is only
/// not shown for partial formulas.
///
bool is_print_combined_result()
{
#ifdef HAVE_MPI
  return !is_print_variables() && is_mpi_master_proc();
#else
  return !is_print_variables();
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

void print_seconds(double seconds)
{
  cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

void print(const string& str)
{
  if (is_print())
    cout << str << endl;
}

void print(const string& str, maxint_t res)
{
  if (is_print())
    cout << str << " = " << res << endl;
}

/// Print result of partial formula
void print(const string& str, maxint_t res, double time)
{
  if (is_print())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << str << " = " << res << endl;
    print_seconds(get_time() - time);
  }
}

/// Used by pi_lmo(x), pi_deleglise_rivat(x)
void print(maxint_t x, int64_t y, int64_t z, int64_t c, int threads)
{
  if (is_print())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    print_threads(threads);
  }
}

/// Only enabled for partial formulas
void print_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
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

/// Only enabled for partial formulas
void print_vars(maxint_t x, int64_t y, int64_t c, int threads)
{
  if (is_print_variables())
  {
    int64_t z = (int64_t)(x / y);
    print(x, y, z, c, threads);
    cout << endl;
  }
}

/// Used by pi_gourdon(x)
void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads)
{
  if (is_print())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "k = " << k << endl;
    cout << "x_star = " << get_x_star_gourdon(x, y) << endl;
    cout << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    cout << "alpha_z = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    print_threads(threads);
  }
}

/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "alpha_y = " << fixed << setprecision(3) << get_alpha_y(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }
}

/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int64_t z, int64_t k, int threads)
{
  if (is_print_variables())
  {
    print_gourdon_vars(x, y, z, k, threads);
    cout << endl;
  }
}

} // namespace
