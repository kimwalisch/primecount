///
/// @file  util.cpp
///        This file contains helper functions and global variables
///        that are initialized with default settings.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <FactorTable.hpp>
#include <int128_t.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <stdint.h>

namespace {

int status_precision_ = -1;

/// Recompute pi(x) with an alternative alpha tuning factor to
/// verify the first result. This redundancy helps guard
/// against potential bugs in primecount: if an error exists,
/// it is highly unlikely that both pi(x) computations would
/// produce the same (incorrect) result.
///
bool double_check_ = false;

/// Truncate a floating point number to 3 digits after the decimal
/// point. This function is used to limit the number of digits after
/// the decimal point of the alpha tuning factor in order to make
/// it more convenient for the user to e.g. type the alpha tuning
/// factor as a command-line parameter.
///
double truncate3(double n)
{
  ASSERT(n <= 1e15);
  return int64_t(n * 1000) / 1000.0;
}

} // namespace

namespace primecount {

int get_status_precision(maxint_t x, int threads)
{
  // Use default precision when no command-line precision provided
  if (status_precision_ < 0 &&
      double(x) > 3e19)
  {
    // D(3e19) takes about 100 seconds using a single thread.
    // Multiplying x by 10 quadruples the runtime, while
    // the speedup is proportional to the number of threads.
    double seconds = 100 * std::pow(double(x) / 3e19, std::log10(4));
    seconds /= max(threads, 1);

    if (seconds >= 10000)
      return 2;
    if (seconds >= 100)
      return 1;
  }

  return max(status_precision_, 0);
}

void set_status_precision(int precision)
{
  status_precision_ = in_between(0, precision, 5);
}

/// Get the time in seconds.
/// Note that according to the documentation of
/// std::chrono::steady_clock: "This clock is not related to wall
/// clock time (for example, it can be time since last reboot)".
/// Hence time will always be fairly small and there won't be any
/// precision issues if we convert the time to double.
///
double get_time()
{
  auto time = std::chrono::steady_clock::now().time_since_epoch();
  auto seconds = std::chrono::duration<double>(time).count();
  // Ensure we have at least microsecond accuracy
  ASSERT(seconds < (1ll << 52) / 1e6);
  return seconds;
}

void set_double_check(bool enable)
{
  double_check_ = enable;
}

/// Tuning factor used in the Lagarias-Miller-Odlyzko,
/// Deleglise-Rivat and Gourdon algorithms.
///
double get_alpha(maxint_t x, int64_t y)
{
  // y = x13 * alpha, thus alpha = y / x13
  double x13 = (double) iroot<3>(x);
  double alpha = (double) y / x13;
  double max_alpha = (double) y;
  int64_t verify_y = (int64_t)(x13 * alpha);

  // Prevent x^(1/3) * alpha = 23.99999...
  if (verify_y < y)
    alpha = std::nextafter(alpha, max_alpha);

  return alpha;
}

/// Get the Lagarias-Miller-Odlyzko alpha tuning factor.
/// alpha = a log(x)^2 + b log(x) + c
/// a, b and c have been determined empirically.
/// @see doc/alpha-factor-lmo.pdf
///
double get_alpha_lmo(maxint_t x, double alpha)
{
  double x16 = (double) iroot<6>(x);

  // use default alpha if no command-line alpha provided
  if (alpha < 1)
  {
    double a = 0.001103;
    double b = -0.00896211;
    double c = 1.00404;
    double logx = std::log((double) x);
    double logx2 = logx * logx;
    alpha = a * logx2 + b * logx + c;
  }

  // Recompute pi(x) with an alternative alpha tuning
  // factor to verify the first result.
  if (double_check_)
    alpha *= 0.97;

  // Preserve 3 digits after decimal point
  alpha = in_between(1, alpha, x16);
  alpha = truncate3(alpha);

  return in_between(1, alpha, x16);
}

/// Get the Deleglise-Rivat alpha tuning factor.
/// alpha = a log(x)^3 + b log(x)^2 + c log(x) + d
/// a, b, c and d have been determined empirically.
/// @see doc/alpha-factor-dr.pdf
///
double get_alpha_deleglise_rivat(maxint_t x, double alpha)
{
  double x16 = (double) iroot<6>(x);

  // Use default alpha
  if (alpha < 1)
  {
    // For x <= 10^9 our default formula does not
    // generate good alpha values. Hence we use
    // another formula optimized for small values.
    if (x <= 1e9)
    {
      double a = 0.078173;
      double b = 1;
      double logx = std::log((double) x);
      alpha = a * logx + b;
    }
    else
    {
      double a = 0.00148918;
      double b = -0.0691909;
      double c = 1.00165;
      double d = 0.372253;
      double logx = std::log((double) x);
      double logx2 = logx * logx;
      double logx3 = logx * logx * logx;
      alpha = a * logx3 + b * logx2 + c * logx + d;
    }
  }

  // Recompute pi(x) with an alternative alpha tuning
  // factor to verify the first result.
  if (double_check_)
    alpha *= 0.97;

  // Preserve 3 digits after decimal point
  alpha = in_between(1, alpha, x16);
  alpha = truncate3(alpha);

  return in_between(1, alpha, x16);
}

/// Get the Xavier Gourdon alpha tuning factor.
/// alpha = a log(x)^3 + b log(x)^2 + c log(x) + d
/// a, b, c and d have been determined empirically.
/// @see doc/alpha-factor-gourdon.pdf
///
double get_alpha_gourdon(maxint_t x, double alpha)
{
  double x16 = (double) iroot<6>(x);

  // Use default alpha
  if (alpha < 1)
  {
    // For x <= 3e11 our default formula does not
    // generate good alpha values. Hence we use
    // another formula optimized for small values.
    if (x <= 3e11)
    {
      double a = 0.156346;
      double b = 1;
      double logx = std::log((double) x);
      alpha = a * logx + b;
    }
    else
    {
      double a = 0.021078652642134526;
      double b = -2.0334362096119443;
      double c = 67.97860096709101;
      double d = -758.0125793471067;
      double logx = std::log((double) x);
      double logx2 = logx * logx;
      double logx3 = logx * logx * logx;
      alpha = a * logx3 + b * logx2 + c * logx + d;

      // If y = x^(1/3) * alpha is too large for 16-bit factor
      // table entries, D uses 32-bit factor table entries
      // which doubles its memory usage. Hence we reduce alpha
      // by 30%, according to our benchmarks this causes a
      // slowdown of at most 3%.
      double x13 = (double) iroot<3>(x);
      int64_t y = (int64_t) (x13 * alpha);
      if (y > FactorTable<uint16_t>::max())
        alpha *= 0.7;
    }
  }

  // Recompute pi(x) with an alternative alpha tuning
  // factor to verify the first result.
  if (double_check_)
    alpha *= 0.97;

  // Preserve 3 digits after decimal point
  alpha = in_between(1, alpha, x16);
  alpha = truncate3(alpha);

  return in_between(1, alpha, x16);
}

/// x_star = max(x^(1/4), x / y^2)
///
/// Gourdon's formulas require x_star <= (x / y)^(1/2):
/// Sigma5 and parts of Sigma2 sum over the primes
/// (x / y)^(1/2) < p <= x^(1/3) independently of x_star, whereas
/// the A, C, Sigma3, Sigma4, Sigma6 formulas and the other
/// parts of Sigma2 depend on x_star. For x^(1/3) < y <= x^(1/2)
/// this bound holds if x / y^2 is rounded down.
///
int64_t get_x_star_gourdon(maxint_t x, int64_t y)
{
  y = max(y, 1);
  maxint_t yy = (maxint_t) y * y;
  int64_t x14 = iroot<4>(x);
  int64_t x_star = max(x / yy, x14);

  // x_star <= y holds if y >= x^(1/3), but for x = 2 and 3
  // we use y = 1. Without this bound x_star > y would cause
  // out of bounds PiTable accesses in AC and Sigma.
  x_star = min(x_star, y);
  x_star = max(x_star, 1);

  return x_star;
}

/// Quickly verify a pi(x) result.
/// Note that this check can only detect miscalculations if the
/// pi(x) result is off by >= sqrt(x) * log(x) / 8π.
///
/// Since we have an extensive test suite that likely finds most
/// implementation bugs, we expect this verification check to
/// mainly detect miscalculations due to hardware errors, such as
/// malfunctioning RAM sticks or PC overclocking issues.
///
void verify_pix(string_view_t pix_function,
                maxint_t x,
                maxint_t pix)
{
  if (x < 2657)
    return;

  maxint_t Lix = Li(x);
  double logx = std::log(x);
  double sqrtx = std::sqrt(x);
  constexpr double PI = 3.14159265358979323846;

  // Lowell Schoenfeld, "Sharper bounds for the Chebyshev functions
  // 𝜃(𝑥) and 𝜓(𝑥). II", Math. Comp., v. 30, 1976.
  if (std::abs(double(pix - Lix)) >= (sqrtx * logx) / (8 * PI))
  {
    std::ostringstream oss;
    oss << "\rprimecount error: " << pix_function << "(" << x << ") = "  << pix  << std::endl
        << "Li(x) = " << Lix << ", sqrt(x) = " << sqrtx << ", log(x) = " << logx << std::endl
        << "Assertion failed: |pi(x) - Li(x)| < sqrt(x) * log(x) / (8 * PI)\n"   << std::endl;
    std::cerr << oss.str();
    std::abort();
  }
}

#if defined(ENABLE_ASSERT)

/// Custom assertion failure handler.
/// On MinGW GCC 16, the standard assert() implementation no longer
/// guarantees a non-returning failure path, which can cause false
/// positive compiler warnings such as -Warray-bounds in debug builds.
/// https://github.com/mingw-w64/mingw-w64/commit/ecf2328a328d11dec7044b40b2b5e93b5b2b9d9e
///
/// Using our own [[noreturn]] handler lets the compiler correctly
/// infer that execution cannot continue after a failed assertion.
///
[[noreturn]]
void assert_failed(const char* assertion,
                   const char* file,
                   unsigned int line,
                   const char* function)
{
  std::string msg("\n");
  msg += std::string(file) + ":" + std::to_string(line);
  msg += ": " + std::string(function);
  msg += ": Assertion `";
  msg += assertion + std::string("' failed.\n\n");

  std::cerr << msg;

  std::abort();
}

#endif

} // namespace
