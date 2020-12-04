///
/// @file  Status.cpp
/// @brief The Status class is used to print the status (in percent)
///        of the formulas related to special leaves. It is used by
///        the D, S2_easy and S2_hard formulas.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Status.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <iostream>
#include <algorithm>
#include <iomanip>

#if defined(_OPENMP)
  #include <omp.h>
#endif

using namespace std;
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

Status::Status(maxint_t x)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

bool Status::isPrint(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

void Status::print(double percent)
{
  double old = percent_;

  if ((percent - old) >= epsilon_)
  {
    percent_ = percent;
    cout << "\rStatus: " << fixed << setprecision(precision_)
         << percent << "%" << flush;
  }
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancer.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
double Status::getPercent(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx)
{
  double p1 = skewed_percent(sum, sum_approx);
  double p2 = skewed_percent(low, limit);

  // When p2 is larger then p1 it is
  // always much more accurate.
  if (p2 > p1)
    return p2;

  // Below 20% p1 is better
  // Above 70% p2 is better
  double c1 = 150 / max(p1, 1.0);
  c1 = in_between(10, c1, 4);
  double c2 = 10 - c1;
  double percent = (c1*p1 + c2*p2) / 10;

  return percent;
}

/// This method is used by S2_hard() and D().
/// This method does not use a lock to synchronize threads
/// as it is only used inside of a critical section inside
/// LoadBalancer.cpp and hence it can never be accessed
/// simultaneously from multiple threads.
///
void Status::print(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx)
{
  double time = get_time();

  if (isPrint(time))
  {
    time_ = time;
    double percent = getPercent(low, limit, sum, sum_approx);
    print(percent);
  }
}

/// Used by S2_easy
void Status::print(int64_t b, int64_t max_b)
{
  // check --status option used
  if (!is_print())
    return;

#if defined(_OPENMP)
  // In order to prevent data races only one thread at a time
  // can enter this code section. In order to make sure that
  // our code scales well up to a very large number of CPU
  // cores, we don't want to use any thread synchronization!
  // In order to achieve this only one of the threads (the
  // main thread) is allowed to print the status, while all
  // other threads do nothing.
  if (omp_get_thread_num() != 0)
    return;
#endif

  double time = get_time();

  if (isPrint(time))
  {
    time_ = time;
    double percent = skewed_percent(b, max_b);
    print(percent);
  }
}

} // namespace
