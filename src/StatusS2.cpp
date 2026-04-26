///
/// @file  StatusS2.cpp
/// @brief The StatusS2 class is used to print the status (in percent)
///        of the formulas related to special leaves. It is used by
///        the D, S2_easy and S2_hard formulas.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusS2.hpp>
#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <algorithm>
#include <cmath>
#include <string>

using namespace primecount;

namespace {

double log_percent(double r, double factor)
{
  double percent = 100.0 * std::log1p(factor * r) / std::log1p(factor);
  return in_between(0.0, percent, 100.0);
}

double blend(double a, double b, double weight_b)
{
  double percent = a * (1.0 - weight_b) + b * weight_b;
  return in_between(0.0, percent, 100.0);
}

double smoothstep(double x)
{
  x = in_between(0.0, x, 1.0);
  return x * x * (3.0 - 2.0 * x);
}

double log_percent(double r,
                   double early_factor,
                   double base_factor,
                   double delay,
                   double cap,
                   double cutoff)
{
  double base = log_percent(r, base_factor);
  double boost = log_percent(r, early_factor);
  boost -= delay * (1.0 - smoothstep(r / cutoff));
  boost = in_between(0.0, boost, cap);

  double percent = std::max(base, boost);
  double floor = std::min(500.0 * r, 0.5);
  floor = in_between(0.0, floor, 100.0);

  return std::max(percent, floor);
}

} // namespace

namespace primecount {

StatusS2::StatusS2(maxint_t x,
                   int64_t y,
                   bool is_print)
{
  double log_y = std::log(y);
  double log10_x = std::log10(double(x));
  y_log_y_ = int64_t(y * log_y);
  x_tune_ = in_between(0.0, (log10_x - 20.0) / 2.0, 1.0);

  if (is_print)
  {
    epsilon_ = 1.0;
    precision_ = get_status_precision(x);

    for (int i = 0; i < precision_; i++)
      epsilon_ /= 10.0;
  }
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
double StatusS2::getPercent(int64_t low, int64_t limit) const
{
  // Works best for >= 90%
  double percent1 = get_percent(low, limit);

  // Works well for <= 20%
  int64_t limit2 = (y_log_y_ * 2) / 3;
  double percent2 = get_percent(low, limit2);
  percent2 = std::min(percent2, 20.0);

  // Works well for >= 20%
  double r = percent1 / 100;
  double small = log_percent(r, 2643.010656, 21.015052, 11.846115, 56.508811, 0.000203036);
  double large = log_percent(r, 25589.45108, 15.357592, 42.898382, 54.704957, 0.000411627);
  double percent3 = blend(small, large, x_tune_);

  double percent23 = std::max(percent2, percent3);

  if (percent23 < 90)
    return percent23;
  else
    return std::min(percent1, percent23);
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
double StatusS2::getStatus(int64_t low, int64_t limit)
{
  double time = get_time();
  double old = time_;

  if ((time - old) >= threshold_)
  {
    time_ = time;
    double old = percent_;
    double percent = getPercent(low, limit);
    if ((percent - old) >= epsilon_)
    {
      percent_ = percent;
      return percent;
    }
  }

  return -1;
}

void StatusS2::print(double percent) const
{
  std::string status = "Status: ";
  status += to_string(percent, precision_);
  status += '%';
  print_status(status);
}

/// This method is used by S2_easy().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// S2_easy.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
void StatusS2::print(int64_t b, int64_t max_b)
{
  double time = get_time();
  double old = time_;

  if ((time - old) >= threshold_)
  {
    time_ = time;

    double p1 = get_percent(b, max_b);
    double p2 = p1 * p1;
    double p3 = p1 * p2;
    double p4 = p2 * p2;
    double c1 = 3.70559815037356886459;
    double c2 = 0.07330455122609925077;
    double c3 = 0.00067895345810494585;
    double c4 = 0.00000216467760881310;
    double percent = -c4*p4 + c3*p3 - c2*p2 + c1*p1;
    percent = in_between(0, percent, 100);

    if ((percent - percent_) >= epsilon_)
    {
      percent_ = percent;
      print(percent);
    }
  }
}

} // namespace
