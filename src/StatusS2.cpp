///
/// @file  StatusS2.cpp
/// @brief The StatusS2 class is used to print the status (in percent)
///        of the formulas related to special leaves. It is used by
///        the D, S2_easy and S2_hard formulas.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusS2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace primecount;

namespace {

/// Since the distribution of the special leaves is highly skewed
/// we cannot simply calculate the percentage of the current
/// computation using the standard linear formula. Hence we use
/// a polynomial formula that grows faster when the value is
/// small and slower towards the end (100%).
/// @see scripts/status_curve_fitting.cpp
///
template <typename T>
double skewed_percent(T x, T y)
{
  double p1 = get_percent(x, y);
  double p2 = p1 * p1;
  double p3 = p1 * p2;
  double p4 = p2 * p2;

  double c1 = 3.70559815037356886459;
  double c2 = 0.07330455122609925077;
  double c3 = 0.00067895345810494585;
  double c4 = 0.00000216467760881310;

  double percent = -c4*p4 + c3*p3 - c2*p2 + c1*p1;
  percent = in_between(0, percent, 100);
  return percent;
}

} // namespace

namespace primecount {

StatusS2::StatusS2(maxint_t x)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
double StatusS2::getPercent(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx)
{
  double p1 = skewed_percent(sum, sum_approx);
  double p2 = skewed_percent(low, limit);

  // When p2 is larger then p1 it is
  // always much more accurate.
  if (p2 > p1)
    return p2;

  // Below 20% p1 is better
  // Above 70% p2 is better
  double c1 = 150 / std::max(p1, 1.0);
  c1 = in_between(10, c1, 4);
  double c2 = 10 - c1;
  double percent = (c1*p1 + c2*p2) / 10;

  return percent;
}

void StatusS2::print(double percent)
{
  double old = percent_;

  if ((percent - old) >= epsilon_)
  {
    percent_ = percent;
    std::cout << "\rStatus: " << std::fixed << std::setprecision(precision_)
              << percent << "%" << std::flush;
  }
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
void StatusS2::print(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx)
{
  double time = get_time();
  double old = time_;

  if ((time - old) >= threshold_)
  {
    time_ = time;
    double percent = getPercent(low, limit, sum, sum_approx);
    print(percent);
  }
}

/// Used by S2_easy.
/// The calling code has to ensure that only 1 thread at a
/// time executes this method.
///
void StatusS2::print(int64_t b, int64_t max_b)
{
  double time = get_time();
  double old = time_;

  if ((time - old) >= threshold_)
  {
    time_ = time;
    double percent = skewed_percent(b, max_b);
    print(percent);
  }
}

} // namespace
