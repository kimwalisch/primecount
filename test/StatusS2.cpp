///
/// @file   StatusS2.cpp
/// @brief  Test the work model used for printing the status of
///         the D, S2_hard and LMO formulas, and test that the
///         sieve_factor of each algorithm is 1 when using the
///         default alpha tuning factor.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusS2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

/// The percentage must be 0 at the start, 100 at the end
/// and monotonically increase in between.
bool is_valid(const StatusS2& status, int64_t limit, double leaves_low)
{
  if (status.getPercent(0) != 0 ||
      std::abs(status.getPercent(limit) - 100) > 1e-9 ||
      status.getPercent(limit + limit / 2) != 100)
    return false;

  double old[3] = { 0, 0, 0 };

  for (int i = 0; i <= 1000; i++)
  {
    // Linear, logarithmic and near the start of the leaves
    double positions[3] = { limit * (i / 1000.0),
                            std::pow(double(limit), i / 1000.0),
                            leaves_low * (1 + i / 200.0) };

    for (int j = 0; j < 3; j++)
    {
      double pos = std::min(positions[j], double(limit));
      double percent = status.getPercent(int64_t(pos));
      if (!std::isfinite(percent) ||
          percent < 0 ||
          percent > 100 ||
          percent < old[j])
        return false;
      old[j] = percent;
    }
  }

  return true;
}

int main()
{
  int64_t limits[3] = { 1000, 46415773936ll, 1ll << 62 };
  double sieve_factors[4] = { 0.001, 1, 50, 1000 };

  for (int64_t limit : limits)
  {
    for (double leaves_low : { 1.0, std::sqrt(double(limit)), limit / 2.0, double(limit) })
    {
      for (double sieve_factor : sieve_factors)
      {
        StatusS2 status(limit, leaves_low, limit, sieve_factor, false);
        std::cout << "StatusS2(limit = " << limit << ", leaves_low = " << leaves_low
                  << ", sieve_factor = " << sieve_factor << ")";
        check(is_valid(status, limit, leaves_low));
      }
    }
  }

  // A larger sieve_factor moves the percentage
  // towards the linear percentage of the sieving.
  int64_t limit = 46415773936ll;
  double leaves_low = 2511886;
  StatusS2 status1(limit, leaves_low, limit, 1, false);
  StatusS2 status2(limit, leaves_low, limit, 50, false);

  for (int i = 1; i < 100; i++)
  {
    int64_t low = int64_t(std::pow(double(limit), i / 100.0));
    double linear = 100.0 * low / limit;
    double diff1 = std::abs(status1.getPercent(low) - linear);
    double diff2 = std::abs(status2.getPercent(low) - linear);
    std::cout << "getPercent(" << low << ") closer to linear for larger sieve_factor";
    check(diff2 <= diff1 + 1e-9);
  }

  // The sieve_factor is computed using get_alpha_*(x) / get_alpha(x, y)
  // which must be 1 when using the default alpha tuning factor.
  for (int i = 12; i <= 18; i++)
  {
    int64_t x = 1;
    for (int j = 0; j < i; j++)
      x *= 10;

    double x13 = (double) iroot<3>(x);

    double alpha = get_alpha_gourdon(x);
    int64_t y = (int64_t) (x13 * alpha);
    y = std::max(y, iroot<3>(x) + 1);
    y = std::min(y, isqrt(x));
    std::cout << "Gourdon: get_alpha_gourdon(" << x << ") / get_alpha(x, y) = " << alpha / get_alpha(x, y);
    check(std::abs(alpha / get_alpha(x, y) - 1) < 1e-3);

    alpha = get_alpha_deleglise_rivat(x);
    y = (int64_t) (x13 * alpha);
    std::cout << "Deleglise-Rivat: get_alpha_deleglise_rivat(" << x << ") / get_alpha(x, y) = " << alpha / get_alpha(x, y);
    check(std::abs(alpha / get_alpha(x, y) - 1) < 1e-3);

    alpha = get_alpha_lmo(x);
    y = (int64_t) (x13 * alpha);
    std::cout << "LMO: get_alpha_lmo(" << x << ") / get_alpha(x, y) = " << alpha / get_alpha(x, y);
    check(std::abs(alpha / get_alpha(x, y) - 1) < 1e-3);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
