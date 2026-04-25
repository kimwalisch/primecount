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
#include <cstdlib>
#include <string>

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

double curve01(double percent)
{
  return in_between(0.0, percent, 100.0);
}

double linear_ratio(int64_t low, int64_t limit)
{
  double percent = get_percent(low, limit);
  return percent / 100.0;
}

double power_percent(double r, double exponent)
{
  return curve01(100.0 * std::pow(r, exponent));
}

double log_percent(double r, double factor)
{
  return curve01(100.0 * std::log1p(factor * r) / std::log1p(factor));
}

double blend(double a, double b, double weight_b)
{
  return curve01(a * (1.0 - weight_b) + b * weight_b);
}

double smoothstep(double x)
{
  x = in_between(0.0, x, 1.0);
  return x * x * (3.0 - 2.0 * x);
}

double delayed(double percent, double r, double delay, double cutoff)
{
  double w = smoothstep(r / cutoff);
  return curve01(percent - delay * (1.0 - w));
}

double early_floor(double percent, double r, double scale, double cap)
{
  double floor = std::min(scale * r, cap);
  return std::max(percent, curve01(floor));
}

double percent_power(double percent, double exponent)
{
  double r = curve01(percent) / 100.0;
  return power_percent(r, exponent);
}

double capped_log_boost_percent(double r,
                                double early_factor,
                                double base_factor,
                                double delay,
                                double cap,
                                double cutoff)
{
  double base = log_percent(r, base_factor);
  double boost = log_percent(r, early_factor);
  boost -= delay * (1.0 - smoothstep(r / cutoff));
  boost = curve01(std::min(boost, cap));
  double percent = std::max(base, boost);

  return early_floor(percent, r, 500.0, 0.5);
}

double tuned_log_percent(double r,
                         double x_tune,
                         double small_early_factor,
                         double small_base_factor,
                         double small_delay,
                         double small_cap,
                         double small_cutoff,
                         double large_early_factor,
                         double large_base_factor,
                         double large_delay,
                         double large_cap,
                         double large_cutoff)
{
  double small = capped_log_boost_percent(r,
                                          small_early_factor,
                                          small_base_factor,
                                          small_delay,
                                          small_cap,
                                          small_cutoff);
  double large = capped_log_boost_percent(r,
                                          large_early_factor,
                                          large_base_factor,
                                          large_delay,
                                          large_cap,
                                          large_cutoff);

  return blend(small, large, x_tune);
}

/// Select one of several experimental no-sum progress curves.
/// This is intentionally controlled through an environment variable
/// so different candidates can be tested without rebuilding.
int get_status_s2_method()
{
  static int method = []()
  {
    const char* env = std::getenv("PRIMECOUNT_S2_STATUS_METHOD");
    if (!env)
      return 0;

    return in_between(0, std::atoi(env), 14);
  }();

  return method;
}

} // namespace

namespace primecount {

StatusS2::StatusS2(maxint_t x)
{
  precision_ = get_status_precision(x);
  x_tune_ = std::log10((double) std::max(x, (maxint_t) 1));
  x_tune_ = in_between(0.0, (x_tune_ - 20.0) / 2.0, 1.0);
  epsilon_ = 1.0;
  for (int i = 0; i < precision_; i++)
    epsilon_ /= 10.0;
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
double StatusS2::getPercent(int64_t low, int64_t limit) const
{
  double r = linear_ratio(low, limit);
  double legacy = skewed_percent(low, limit);

  switch (get_status_s2_method())
  {
    case  0: return tuned_log_percent(r, x_tune_,
                                      2643.010656, 21.015052, 11.846115, 56.508811, 0.000203036,
                                      25589.451080, 15.357592, 42.898382, 54.704957, 0.000411627);
    case  1: return tuned_log_percent(r, x_tune_,
                                      3071.222874, 18.992829, 18.456959, 57.874339, 0.000425254,
                                      30338.257157, 15.752385, 44.494719, 52.798611, 0.000413516);
    case  2: return tuned_log_percent(r, x_tune_,
                                      2171.215935, 14.132140, 10.683080, 52.926555, 0.001359935,
                                      23855.906060, 13.437586, 53.480088, 50.950426, 0.000395379);
    case  3: return capped_log_boost_percent(r, 2643.010656, 21.015052, 11.846115, 56.508811, 0.000203036);
    case  4: return capped_log_boost_percent(r, 25589.451080, 15.357592, 42.898382, 54.704957, 0.000411627);
    case  5: return early_floor(delayed(log_percent(r, 2000.0), r, 8.0, 0.006), r, 500.0, 0.5);
    case  6: return early_floor(delayed(log_percent(r, 2000.0), r, 6.0, 0.006), r, 300.0, 0.5);
    case  7: return early_floor(delayed(log_percent(r, 2000.0), r, 7.0, 0.006), r, 300.0, 0.5);
    case  8: return early_floor(delayed(log_percent(r, 2000.0), r, 8.0, 0.006), r, 500.0, 0.5);
    case  9: return log_percent(r, 500.0);
    case 10: return log_percent(r, 1000.0);
    case 11: return log_percent(r, 2000.0);
    case 12: return percent_power(legacy, 0.35);
    case 13: return capped_log_boost_percent(r, 3071.222874, 18.992829, 18.456959, 57.874339, 0.000425254);
    case 14: return capped_log_boost_percent(r, 30338.257157, 15.752385, 44.494719, 52.798611, 0.000413516);
    default:
      return legacy;
  }
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
    double old = percent_;
    double percent = skewed_percent(b, max_b);
    if ((percent - old) >= epsilon_)
    {
      percent_ = percent;
      print(percent);
    }
  }
}

} // namespace
