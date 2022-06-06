///
/// @file  util.cpp
///        This file contains helper functions and global variables
///        that are initialized with default settings.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <calculator.hpp>
#include <int128_t.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <to_string.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <stdint.h>
#include <utility>

using std::min;
using std::max;

namespace {

int status_precision_ = -1;

// Tuning factor used in the Lagarias-Miller-Odlyzko
// and Deleglise-Rivat algorithms.
double alpha_ = -1;

// Tuning factor used in Xavier Gourdon's algorithm
double alpha_y_ = -1;

// Tuning factor used in Xavier Gourdon's algorithm
double alpha_z_ = -1;

/// Truncate a floating point number to 3 digits after the decimal
/// point. This function is used limit the number of digits after the
/// decimal point of the alpha tuning factor in order to make it more
/// convenient for the user to e.g. type the alpha tuning factor as
/// a command-line parameter.
///
double truncate3(double n)
{
  return (int64_t)(n * 1000) / 1000.0;
}

} // namespace

namespace primecount {

/// The compiler supports the non standard __int128_t type, but the
/// standard int128_t type is missing in <stdint.h>. We need to
/// define a few functions that are not supported by the C++ STL.
///
#if defined(HAVE_NON_STANDARD__INT128_T)

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
    return "-" + to_string((uint128_t) -n);
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

maxint_t to_maxint(const std::string& expr)
{
  // Check if n <= max, if n is only
  // composed of regular digits.
  if (expr.find_first_not_of("0123456789") == std::string::npos)
  {
    // Remove leading zeros
    std::size_t pos = expr.find_first_not_of("0");

    if (pos != std::string::npos)
    {
      std::string n = expr.substr(pos);
      maxint_t limit = std::numeric_limits<maxint_t>::max();
      std::string max_n = to_string(limit);

      if (n.size() > max_n.size() ||
        (n.size() == max_n.size() && n > max_n))
      {
        std::string msg = "number too large: " + n;
        throw primecount_error(msg);
      }
    }
  }

  return calculator::eval<maxint_t>(expr);
}

int ideal_num_threads(int threads, int64_t sieve_limit, int64_t thread_threshold)
{
  thread_threshold = max((int64_t) 1, thread_threshold);
  int64_t max_threads = ceil_div(sieve_limit, thread_threshold);
  return in_between(1, threads, max_threads);
}

int get_status_precision(maxint_t x)
{
  // use default precision when no command-line precision provided
  if (status_precision_ < 0)
  {
    if ((double) x >= 1e23)
      return 2;
    if ((double) x >= 1e21)
      return 1;
  }

  return max(status_precision_, 0);
}

void set_status_precision(int precision)
{
  status_precision_ = in_between(0, precision, 5);
}

/// Get the time in seconds (with microsecond accuracy).
/// Note that according to the documentation of
/// std::chrono::steady_clock: "This clock is not related to wall
/// clock time (for example, it can be time since last reboot)".
/// Hence time will always be fairly small and there won't be any
/// precision issues if we convert the time to double.
///
double get_time()
{
  auto now = std::chrono::steady_clock::now();
  auto time = now.time_since_epoch();
  auto micro = std::chrono::duration_cast<std::chrono::microseconds>(time);
  ASSERT(micro.count() < (1ll << 52));
  return (double) micro.count() / 1e6;
}

void set_alpha(double alpha)
{
  // If alpha < 1 then we compute a good
  // alpha tuning factor at runtime.
  if (alpha < 1.0)
    alpha_ = -1;
  else
    alpha_ = truncate3(alpha);
}

void set_alpha_y(double alpha_y)
{
  // If alpha_y < 1 then we compute a good
  // alpha tuning factor at runtime.
  if (alpha_y < 1.0)
    alpha_y_ = -1;
  else
    alpha_y_ = truncate3(alpha_y);
}

void set_alpha_z(double alpha_z)
{
  // If alpha_z < 1 then we compute a good
  // alpha tuning factor at runtime.
  if (alpha_z < 1.0)
    alpha_z_ = -1;
  else
    alpha_z_ = truncate3(alpha_z);
}

/// Tuning factor used in the Lagarias-Miller-Odlyzko
/// and Deleglise-Rivat algorithms.
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

/// Tuning factor used in Xavier Gourdon's algorithm.
double get_alpha_y(maxint_t x, int64_t y)
{
  // y = x13 * alpha_y, thus alpha = y / x13
  double x13 = (double) iroot<3>(x);
  double alpha_y = (double) y / x13;
  double max_alpha_y = (double) y;
  int64_t verify_y = (int64_t)(x13 * alpha_y);

  // Prevent x^(1/3) * alpha_y = 23.99999...
  if (verify_y < y)
    alpha_y = std::nextafter(alpha_y, max_alpha_y);

  return alpha_y;
}

/// Tuning factor used in Xavier Gourdon's algorithm.
double get_alpha_z(int64_t y, int64_t z)
{
  // z = y * alpha_z, thus alpha_z = z / y
  double alpha_z = (double) z / y;
  double max_alpha_z = (double) z;
  int64_t verify_z = (int64_t)(y * alpha_z);

  // Prevent y * alpha_z = 23.99999...
  if (verify_z < z)
    alpha_z = std::nextafter(alpha_z, max_alpha_z);

  return alpha_z;
}

/// Get the Lagarias-Miller-Odlyzko alpha tuning factor.
/// alpha = a log(x)^2 + b log(x) + c
/// a, b and c have been determined empirically.
/// @see doc/alpha-factor-lmo.pdf
///
double get_alpha_lmo(maxint_t x)
{
  double alpha = alpha_;
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
double get_alpha_deleglise_rivat(maxint_t x)
{
  double alpha = alpha_;
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

  // Preserve 3 digits after decimal point
  alpha = in_between(1, alpha, x16);
  alpha = truncate3(alpha);

  return in_between(1, alpha, x16);
}

/// In Xavier Gourdon's algorithm there are 2 alpha tuning
/// factors. The alpha_y tuning factor should grow like
/// O(log(x)^3) and the alpha_z tuning factor is a small
/// constant. Both alpha_y and alpha_z should be determined
/// experimentally by running benchmarks.
/// @see doc/alpha-factor-gourdon.pdf
///
/// y = x^(1/3) * alpha_y, with alpha_y >= 1.
/// z = y * alpha_z, with alpha_z >= 1.
/// alpha_y * alpha_z <= x^(1/6)
///
std::pair<double, double> get_alpha_gourdon(maxint_t x)
{
  double alpha_y = alpha_y_;
  double alpha_z = alpha_z_;
  double x16 = (double) iroot<6>(x);
  double logx = std::log((double) x);
  double alpha_yz;

  // For x <= 10^11 our default formula does not
  // generate good alpha values. Hence we use
  // another formula optimized for small values.
  if (x <= 1e11)
  {
    double a = 0.078173;
    double b = 1;
    alpha_yz = a * logx + b;
  }
  else
  {
    double a = 0.00526934;
    double b = -0.495545;
    double c = 16.5791;
    double d = -183.836;
    double logx2 = logx * logx;
    double logx3 = logx * logx * logx;
    alpha_yz = a * logx3 + b * logx2 + c * logx + d;
  }

  // Use default alpha_z
  if (alpha_z < 1)
  {
    // y = x^(1/3) * alpha_y
    // z = y * alpha_z
    //
    // alpha_y should grow like O(log(x)^3) just like in the
    // Deleglise-Rivat algorithm whereas alpha_z is a small tuning
    // factor usually within [1, 4]. In my opinion the algorithm is
    // theoretically most efficient (i.e. uses the fewest number of
    // instructions) if (y == z), hence if alpha_z = 1. Because when
    // setting y to a value smaller than z this will decrease the
    // number of sparse easy leaves (which can be computed more
    // efficiently than other types of leaves) and increase the
    // number of other types of leaves.
    //
    // By setting alpha_z to a value > 1 this will cause y to be set
    // to a value < z which will generally improve the cache
    // efficiency of the algorithm but as a drawback also increase
    // the number of instructions used by the algorithm. The C1
    // algorithm (in AC.cpp) has severe scaling issues above 10^23
    // as it is not segmented and requires frequent thread
    // synchronization. The larger alpha_z, the less work there will
    // be in the C1 algorithm. Hence for computations >= 10^23 using
    // an alpha_z > 1 will likely improve performance.
    alpha_z = 2;

    // alpha_z should be significantly smaller than alpha_y
    alpha_z = in_between(1, alpha_yz / 5, alpha_z);
  }

  // Use default alpha_y
  if (alpha_y < 1)
    alpha_y = alpha_yz / alpha_z;

  // Preserve 3 digits after decimal point
  alpha_y = in_between(1, alpha_y, x16);
  alpha_y = truncate3(alpha_y);
  alpha_z = truncate3(alpha_z);

  // Ensure alpha_y * alpha_z <= x^(1/6)
  alpha_y = in_between(1, alpha_y, x16);
  double max_alpha_z = max(1.0, x16 / alpha_y);
  alpha_z = in_between(1, alpha_z, max_alpha_z);

  return std::make_pair(alpha_y, alpha_z);
}

/// x_star = max(x^(1/4), x / y^2)
///
/// After my implementation of Xavier Gourdon's algorithm worked for
/// the first time there were still many miscalculations mainly for
/// small numbers < 10^6. By debugging I found that most errors were
/// related to the Sigma formulas (Σ0 - Σ6) and the x_star variable
/// was responsible for most errors. For some unknown reason the
/// bounds from Xavier's paper (max(x^(1/4), x / y^2)) don't seem to
/// be enough. By trial and error I figured out a few more bounds that
/// fix all miscalculations in my implementation.
///
int64_t get_x_star_gourdon(maxint_t x, int64_t y)
{
  // For some unknown reason it is necessary
  // to round up (x / y^2). Without rounding up
  // there are many miscalculations below 2000
  // in my implementation.
  y = max(y, (int64_t) 1);
  maxint_t yy = (maxint_t) y * y;
  maxint_t x_div_yy = ceil_div(x, yy);

  int64_t x_star = (int64_t) max(iroot<4>(x), x_div_yy);
  int64_t sqrt_xy = (int64_t) isqrt(x / y);

  // x_star <= y
  // x_star <= (x / y)^(1/2)
  // The bounds above are missing in Xavier Gourdon's
  // paper. Without these bounds many of the 7 Sigma
  // formulas (Σ0 - Σ6) return incorrect results for
  // numbers below 10^6.
  x_star = min(x_star, y);
  x_star = min(x_star, sqrt_xy);
  x_star = max(x_star, (int64_t) 1);

  return x_star;
}

} // namespace
