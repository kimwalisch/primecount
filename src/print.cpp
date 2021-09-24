///
/// @file  print.cpp
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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

namespace {

bool print_ = false;
bool print_variables_ = false;

bool is_print_variables()
{
  return print_variables_;
}

void print_threads(int threads)
{
  std::cout << "threads = " << threads << std::endl;
}

} // naespace

namespace primecount {

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
  return !is_print_variables();
}

void set_print(bool print)
{
  print_ = print;
}

void set_print_variables(bool print_variables)
{
  print_variables_ = print_variables;
}

void print_seconds(double seconds)
{
  std::cout << "Seconds: " << std::fixed << std::setprecision(3) << seconds << std::endl;
}

void print(const std::string& str)
{
  std::cout << str << std::endl;
}

void print(const std::string& str, maxint_t res)
{
  std::cout << str << " = " << res << std::endl;
}

void print(const std::string& str, maxint_t res, double time)
{
  std::cout << "\r" << std::string(50,' ') << "\r";
  std::cout << "Status: 100%" << std::endl;
  std::cout << str << " = " << res << std::endl;
  print_seconds(get_time() - time);
}

/// Used by pi_lmo(x), pi_deleglise_rivat(x)
void print(maxint_t x, int64_t y, int64_t z, int64_t c, int threads)
{
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "z = " << z << std::endl;
  std::cout << "c = " << c << std::endl;
  std::cout << "alpha = " << std::fixed << std::setprecision(3) << get_alpha(x, y) << std::endl;
  print_threads(threads);
}

/// Only enabled for partial formulas
void print_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
  {
    maxint_t z = x / y;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl;
    std::cout << "alpha = " << std::fixed << std::setprecision(3) << get_alpha(x, y) << std::endl;
    print_threads(threads);
    std::cout << std::endl;
  }
}

/// Only enabled for partial formulas
void print_vars(maxint_t x, int64_t y, int64_t c, int threads)
{
  if (is_print_variables())
  {
    int64_t z = (int64_t)(x / y);
    print(x, y, z, c, threads);
    std::cout << std::endl;
  }
}

/// Used by pi_gourdon(x)
void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads)
{
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "z = " << z << std::endl;
  std::cout << "k = " << k << std::endl;
  std::cout << "x_star = " << get_x_star_gourdon(x, y) << std::endl;
  std::cout << "alpha_y = " << std::fixed << std::setprecision(3) << get_alpha_y(x, y) << std::endl;
  std::cout << "alpha_z = " << std::fixed << std::setprecision(3) << get_alpha_z(y, z) << std::endl;
  print_threads(threads);
}

/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
  {
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "alpha_y = " << std::fixed << std::setprecision(3) << get_alpha_y(x, y) << std::endl;
    print_threads(threads);
    std::cout << std::endl;
  }
}

/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int64_t z, int64_t k, int threads)
{
  if (is_print_variables())
  {
    print_gourdon(x, y, z, k, threads);
    std::cout << std::endl;
  }
}

} // namespace
