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
#include <stdint.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;

namespace {

bool print_status_ = false;

bool print_variables_ = false;

}

namespace primecount {

void set_print_status(bool print_status)
{
  print_status_ = print_status;
}

void set_print_variables(bool print_variables)
{
  print_variables_ = print_variables;
}

bool print_result()
{
  return !print_variables();
}

bool print_status()
{
  return print_status_;
}

bool print_variables()
{
  return print_variables_;
}

void print(const string& str)
{
  if (print_status())
    cout << str << endl;

  if (is_log())
  {
    ofstream outfile("primecount.log", std::ofstream::out | std::ofstream::app);

    if (outfile.is_open())
    {
      outfile << str << endl;
      outfile.close();
    }
  }
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

  if (is_log())
  {
    ofstream outfile("primecount.log", std::ofstream::out | std::ofstream::app);

    if (outfile.is_open())
    {
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "c = " << c << endl;
      outfile << "alpha = " << fixed << setprecision(3) << alpha << endl;
      outfile << "threads = " << validate_threads(threads) << endl;

      outfile.close();
    }
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
    cout << "threads = " << validate_threads(threads) << endl;
    cout << endl;
  }

  if (is_log())
  {
    ofstream outfile("primecount.log", std::ofstream::out | std::ofstream::app);

    if (outfile.is_open())
    {
      maxint_t z = x / y;
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
      outfile << "threads = " << validate_threads(threads) << endl;
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
    cout << "threads = " << validate_threads(threads) << endl;
    cout << endl;
  }

  if (is_log())
  {
    ofstream outfile("primecount.log", std::ofstream::out | std::ofstream::app);

    if (outfile.is_open())
    {
      maxint_t z = x / y;
      outfile << "x = " << x << endl;
      outfile << "y = " << y << endl;
      outfile << "z = " << z << endl;
      outfile << "c = " << c << endl;
      outfile << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
      outfile << "threads = " << validate_threads(threads) << endl;
      outfile << endl;

      outfile.close();
    }
  }
}

void print(const string& res_name, maxint_t res, double time)
{
  double seconds = get_wtime() - time;

  if (print_status())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << res_name << " = " << res << endl;
    print_seconds(seconds);
  }

  if (is_log())
  {
    ofstream outfile("primecount.log", std::ofstream::out | std::ofstream::app);

    if (outfile.is_open())
    {
      outfile << "\r" << string(50,' ') << "\r";
      outfile << "Status: 100%" << endl;
      outfile << res_name << " = " << res << endl;
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
