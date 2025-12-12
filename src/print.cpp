///
/// @file  print.cpp
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <stdint.h>

#include <algorithm>
#include <iostream>

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

} // namespace

#if __cplusplus >= 201703L && \
    __has_include(<charconv>)

#include <primecount.hpp>
#include <Vector.hpp>
#include <macros.hpp>

#include <charconv>
#include <cmath>

namespace primecount {

std::string to_string(double x, int precision)
{
  ASSERT(precision >= 0);
  precision = std::min(precision, 10);

  std::chars_format format = (std::abs(x) < 1e16)
    ? std::chars_format::fixed
    : std::chars_format::scientific;

  // The double value 1e16-1 with the maximum
  // precision 10 requires 27 characters.
  Array<char, 32> buffer;

  std::to_chars_result res = std::to_chars(
    buffer.data(),
    buffer.data() + buffer.size(),
    x, format, precision
  );

  if (res.ec == std::errc{})
    return std::string(buffer.data(), res.ptr);
  else
    throw primecount_error("to_string(double, int): conversion failed!");
}

} // namespace

#else

#include <sstream>

namespace primecount {

/// std::ostringstream has poor performance.
/// Hence, we only use it if the compiler does
/// not support std::to_chars() from C++17.
///
std::string to_string(double x, int precision)
{
  std::ostringstream oss;
  oss << std::fixed;
  oss.precision(precision);
  oss << x;

  return oss.str();
}

} // namespace

#endif

namespace primecount {

/// The compiler supports the non standard __int128_t
/// type, but the int128_t type is missing in <stdint.h>.
/// Hence, we need to define a few int128_t functions
/// that are not supported by the C++ STL.
///
#if defined(ENABLE_INT128_TO_STRING)

std::string to_string(uint128_t n)
{
  std::string str;

  while (n > 0)
  {
    str += '0' + n % 10;
    n /= 10;
  }

  if (str.empty())
    str = "0";

  std::reverse(str.begin(), str.end());

  return str;
}

std::string to_string(int128_t n)
{
  if (n >= 0)
    return to_string((uint128_t) n);
  else
  {
    // -n causes undefined behavior for n = INT128_MIN.
    // Hence we use the defined two's complement negation: ~n + 1.
    // Casting ~n to unsigned ensures the result of the addition
    // (2^127 for INT128_MIN) is safely stored in a uint128_t
    // without signed overflow.
    uint128_t abs_n = uint128_t(~n) + 1;
    return "-" + to_string(abs_n);
  }
}

std::ostream& operator<<(std::ostream& stream, int128_t n)
{
  stream << to_string(n);
  return stream;
}

std::ostream& operator<<(std::ostream& stream, uint128_t n)
{
  stream << to_string(n);
  return stream;
}

#endif

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
  std::cout << "Seconds: " << to_string(seconds, 3) << std::endl;
}

void print(string_view_t str)
{
  std::cout << str << std::endl;
}

void print(string_view_t str, maxint_t res)
{
  std::cout << str << " = " << res << std::endl;
}

void print(string_view_t str, maxint_t res, double time)
{
  // We overwrite the current text line,
  // which could be e.g.:
  // "Status: 99.9999999991%"
  // "Segments; 123456789/123456789"
  std::cout << "\rStatus: 100%                                 " << std::endl;
  std::cout << str << " = " << res << std::endl;
  print_seconds(get_time() - time);
}

/// Used by pi_lmo(x), pi_deleglise_rivat(x)
void print(maxint_t x, int64_t y, int64_t z, int64_t c, int threads)
{
  double alpha = get_alpha(x, y);
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "z = " << z << std::endl;
  std::cout << "c = " << c << std::endl;
  std::cout << "alpha = " << to_string(alpha, 3) << std::endl;
  print_threads(threads);
}

/// Only enabled for partial formulas
void print_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
  {
    maxint_t z = x / y;
    double alpha = get_alpha(x, y);
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl;
    std::cout << "alpha = " << to_string(alpha, 3) << std::endl;
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
  int64_t x_star = get_x_star_gourdon(x, y);
  double alpha_y = get_alpha_y(x, y);
  double alpha_z = get_alpha_z(y, z);

  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "z = " << z << std::endl;
  std::cout << "k = " << k << std::endl;
  std::cout << "x_star = " << x_star << std::endl;
  std::cout << "alpha_y = " << to_string(alpha_y, 3) << std::endl;
  std::cout << "alpha_z = " << to_string(alpha_z, 3) << std::endl;

  print_threads(threads);
}

/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int threads)
{
  if (is_print_variables())
  {
    double alpha_y = get_alpha_y(x, y);
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "alpha_y = " << to_string(alpha_y, 3) << std::endl;
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

void print_nth_prime_sieve(uint64_t n,
                           bool sieve_forward,
                           maxint_t nth_prime_approx,
                           uint64_t dist_approx,
                           uint64_t thread_dist,
                           int threads)
{
  std::cout << "n = " << n << std::endl;
  std::cout << "sieve_forward = " << (sieve_forward ? "true" : "false") << std::endl;
  std::cout << "nth_prime_approx = " << nth_prime_approx << std::endl;
  std::cout << "dist_approx = " << dist_approx << std::endl;
  std::cout << "thread_dist = " << thread_dist << std::endl;
  std::cout << "threads = " << threads << std::endl;
}

} // namespace
