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
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#if defined(_OPENMP)
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Since the distribution of the special leaves is highly skewed
/// we cannot simply calculate the percentage of the current
/// computation using the well known linear formula. The
/// implementation below is a hack which skews the percent result
/// in order to get a more accurate estimation of the current
/// computation status.
///
template <typename T>
double skewed_percent(T x, T y)
{
  double exp = 0.96;
  double percent = get_percent(x, y);
  double base = exp + percent / (101 / (1 - exp));
  double low = pow(base, 100.0);
  double dividend = pow(base, percent) - low;
  percent = 100 - (100 * dividend / (1 - low));
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
    ostringstream status;
    status << "\rStatus: " << fixed << setprecision(precision_) << percent << "%";
    cout << status.str() << flush;
  }
}

/// Executed at the beginning of each segment
void StatusAC::init()
{
  if (!is_print())
    return;

#if defined(_OPENMP)
  if (omp_get_thread_num() != 0)
    return;
#endif

  if (percent_segment_ == -1)
  {
    percent_total_ = 0;
    percent_segment_ = 80;
  }
  else
  {
    percent_total_ += percent_segment_;
    percent_segment_ = (100 - percent_total_) / 3;
  }
}

void StatusAC::print(int64_t b, int64_t max_b)
{
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

} // namespace
