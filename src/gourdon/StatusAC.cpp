///
/// @file  StatusAC.cpp
/// @brief The StatusAC class is used to print the status (in percent)
///        of the A & C formulas in Xavier Gourdon's algorithm.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusAC.hpp>
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
  double p5 = p2 * p3;

  double c1 = 5.71586523529789379523;
  double c2 = 0.17456675134934354428;
  double c3 = 0.00283157450325605000;
  double c4 = 0.00002327032697966551;
  double c5 = 0.00000007695391846741;

  double percent = c5*p5 - c4*p4 + c3*p3 - c2*p2 + c1*p1;
  percent = in_between(0, percent, 100);
  return percent;
}

} // namespace

namespace primecount {

StatusAC::StatusAC(maxint_t x)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

bool StatusAC::isPrint(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

void StatusAC::print(double percent)
{
  double old = percent_;

  if ((percent - old) >= epsilon_)
  {
    percent_ = percent;
    cout << "\rStatus: " << fixed << setprecision(precision_)
         << percent << "%" << flush;
  }
}

/// Initialize for next segment
void StatusAC::next()
{
  percent_total_ += percent_segment_;
  percent_segment_ = (100 - percent_total_) / 3.5;
}

void StatusAC::print(int64_t b, int64_t max_b)
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
    percent = percent_total_ + (percent_segment_ / 100 * percent);
    print(percent);
  }
}

double StatusAC::getPercent(int64_t b, int64_t max_b)
{
  double percent = skewed_percent(b, max_b);
  percent = percent_total_ + (percent_segment_ / 100 * percent);
  return percent;
}

} // namespace
