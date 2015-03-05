///
/// @file  print.cpp
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primecount-internal.hpp>
#include <int128.hpp>
#include <pmath.hpp>
#include <stdint.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_status_ = false;

bool partial_computation_ = true;

}

namespace primecount {

void set_print_status(bool print_status)
{
  print_status_ = print_status;
}

void set_partial_computation(bool partial_computation)
{
  partial_computation_ = partial_computation;
}

bool print_result()
{
  return partial_computation();
}

bool print_status()
{
  return print_status_;
}

bool partial_computation()
{
  return partial_computation_;
}

void print(const string& str)
{
  if (print_status())
    cout << str << endl;
}

void print(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
{
  if (print_status())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }
}

void print(maxint_t x, int64_t y, int64_t c, int threads)
{
  if (print_status() && !partial_computation())
  {
    maxint_t z = x / y;
    double alpha = (double) y / (double) iroot<3>(x);

    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    cout << "threads = " << validate_threads(threads) << endl;
    cout << endl;
  }
}

void print(const string& res_name, maxint_t res, double time)
{
  if (print_status())
  {
    cout << "\r" << string(40,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << res_name << " = " << res << endl;
    print_seconds(get_wtime() - time);
  }
}

void print_seconds(double seconds)
{
  if (print_status())
    cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

} // namespace
