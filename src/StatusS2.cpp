///
/// @file  StatusS2.cpp
/// @brief The StatusS2 class is used to print the status (in percent)
///        of the formulas related to special leaves. It is used by
///        the D, S2_easy, S2_hard and LMO formulas.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusS2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <cmath>
#include <string>

namespace primecount {

StatusS2::StatusS2(maxint_t x, int threads, bool is_print)
{
  ASSERT(x >= 1);

  if (is_print)
  {
    epsilon_ = 1.0;
    precision_ = get_status_precision(x, threads);

    for (int i = 0; i < precision_; i++)
      epsilon_ /= 10.0;
  }
}

/// Used by S2_hard(), D() and LMO, which sieve the interval
/// [0, limit[. The work is modeled as evenly distributed
/// sieving work plus the work of the special leaves, which
/// starts near leaves_low and then decreases like O(n^(-3/2)),
/// i.e. a Lomax distribution with shape 1/2. The constants
/// have been determined empirically for 10^16 <= x <= 10^23
/// using the default alpha. sieve_factor adjusts the share of
/// the sieving work for other alphas; D and S2_hard use 1
/// for the default alpha.
///
StatusS2::StatusS2(maxint_t x,
                   double leaves_low,
                   int64_t limit,
                   double sieve_factor,
                   int threads,
                   bool is_print)
  : StatusS2(x, threads, is_print)
{
  ASSERT(sieve_factor > 0);

  limit_ = (double) max(limit, 1);
  leaves_low = in_between(1.0, leaves_low, limit_);
  leaves_low_ = 1.5 * leaves_low;
  leaves_scale_ = 3.7 * leaves_low;
  sieve_work_ = sieve_factor * 0.066 * std::pow(limit_ / leaves_low, 0.16);
  total_work_ = work(limit_);
}

/// Estimated work of the interval [0, low[,
/// the special leaves have a total work of 1.
///
double StatusS2::work(double low) const
{
  double sieve = sieve_work_ * low / limit_;

  if (low <= leaves_low_)
    return sieve;

  double leaves = 1 - 1 / std::sqrt(1 + (low - leaves_low_) / leaves_scale_);
  return sieve + leaves;
}

/// Estimated percentage of the total
/// work of the interval [0, low[.
///
double StatusS2::getPercent(int64_t low) const
{
  double l = double(low);
  double percent = 100 * work(l) / total_work_;
  return in_between(0.0, percent, 100.0);
}

/// This method is used by S2_hard(), D() and LMO.
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancerS2.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
void StatusS2::print_S2_hard(double percent)
{
  double old = percent_;

  if ((percent - old) >= epsilon_)
  {
    percent_ = percent;
    std::string status = "Status: ";
    status += to_string(min(percent, 100 - epsilon_), precision_);
    status += '%';
    print_status(status);
  }
}

/// This method is used by S2_easy().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// S2_easy.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
void StatusS2::print_S2_easy(int64_t b, int64_t max_b)
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
      std::string status = "Status: ";
      status += to_string(percent, precision_);
      status += '%';
      print_status(status);
    }
  }
}

} // namespace
