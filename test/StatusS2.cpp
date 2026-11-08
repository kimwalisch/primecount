///
/// @file   StatusS2.cpp
/// @brief  Test the work model used for printing the status of
///         the D, S2_hard and LMO formulas, and test that the
///         alpha factor of each algorithm's sieve_work is 1 when
///         using the default alpha tuning factor.
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
bool is_valid(const StatusS2& status, int64_t limit, double leaves_start)
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
                            leaves_start * (1 + i / 200.0) };

    for (int j = 0; j < 3; j++)
    {
      double pos = std::min(positions[j], double(limit));
      double percent = status.getPercent(int64_t(pos));
      // Allow rounding differences between x87 registers and stored doubles
      if (!std::isfinite(percent) ||
          percent < 0 ||
          percent > 100 ||
          percent + 1e-9 < old[j])
        return false;
      old[j] = percent;
    }
  }

  return true;
}

int main()
{
  int64_t limits[3] = { 1000, 46415773936ll, 1ll << 62 };
  // The alpha factor (default_alpha / alpha)^exponent of sieve_work
  // ranges from 1.9e-7 (LMO, alpha = x^(1/6)) to 1720 (D, alpha = 1)
  double sieve_factors[6] = { 1e-7, 0.001, 1, 50, 1000, 1e4 };
  WorkModelS2 models[3] = { { 0, 0.0149, 1.9307, 0.41271, 0.34966 },
                            { 0, 0.076243, 4.3158, 0.61257, 0.1403 },
                            { 0, 0.47309, 1.1764, 0.36355, 0.030408 } };

  for (int64_t limit : limits)
  {
    for (double leaves_start : { 1.0, std::sqrt(double(limit)), limit / 2.0, double(limit) })
    {
      for (const WorkModelS2& algo_model : models)
      {
        for (double sieve_factor : sieve_factors)
        {
          WorkModelS2 model = algo_model;
          model.leaves_start = leaves_start;
          model.sieve_work *= sieve_factor;
          StatusS2 status(limit, limit, model, 1, false);
          std::cout << "StatusS2(limit = " << limit << ", leaves_start = " << leaves_start
                    << ", leaves_shape = " << model.leaves_shape
                    << ", sieve_work = " << model.sieve_work << ")";
          check(is_valid(status, limit, leaves_start));
        }
      }
    }
  }

  // A larger sieve_work moves the percentage
  // towards the linear percentage of the sieving.
  int64_t limit = 46415773936ll;

  for (WorkModelS2 model1 : models)
  {
    model1.leaves_start = 3767829;
    WorkModelS2 model2 = model1;
    model2.sieve_work *= 50;
    StatusS2 status1(limit, limit, model1, 1, false);
    StatusS2 status2(limit, limit, model2, 1, false);

    for (int i = 1; i < 100; i++)
    {
      int64_t low = int64_t(std::pow(double(limit), i / 100.0));
      double linear = 100.0 * low / limit;
      double diff1 = std::abs(status1.getPercent(low) - linear);
      double diff2 = std::abs(status2.getPercent(low) - linear);
      std::cout << "getPercent(" << low << ") closer to linear for larger sieve_work";
      check(diff2 <= diff1 + 1e-9);
    }
  }

  // The sieve_work is adjusted using get_alpha_*(x) / get_alpha(x, y)
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
