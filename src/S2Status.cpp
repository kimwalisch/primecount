///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent,
///        requires --status[=N] option.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace primecount {

S2Status::S2Status(maxint_t x)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

double S2Status::getPercent(int64_t low, int64_t limit, maxint_t S2, maxint_t S2_approx)
{
  double p1 = skewed_percent(low, limit);
  double p2 = skewed_percent(S2, S2_approx);
  double percent = max(p1, p2);

  // p2 is just an approximation,
  // use p1 near the end
  if (p2 > 95.0)
    percent = max(p1, 95.0);

  return percent;
}

/// Dirty hack!
double S2Status::skewed_percent(maxint_t x, maxint_t y)
{
  double exp = 0.96;
  double percent = get_percent(x, y);
  double base = exp + percent / (101 / (1 - exp));
  double low = pow(base, 100.0);
  double dividend = pow(base, percent) - low;
  percent = 100 - (100 * dividend / (1 - low));

  return percent;
}

#if defined(_OPENMP)

bool S2Status::is_print(double time)
{
  TryLock lock(lock_);
  if (lock.ownsLock())
  {
    double old = time_;
    return old == 0 ||
          (time - old) >= is_print_;
  }

  return false;
}

#else

bool S2Status::is_print(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

#endif

void S2Status::print(maxint_t n, maxint_t limit)
{
  double time = get_time();

  if (is_print(time))
  {
    time_ = time;

    double percent = skewed_percent(n, limit);
    double old = percent_;

    if ((percent - old) >= epsilon_)
    {
      percent_ = percent;
      ostringstream status;
      ostringstream out;

      status << "Status: " << fixed << setprecision(precision_) << percent << "%";
      size_t spaces = status.str().length();
      string reset_line = "\r" + string(spaces,' ') + "\r";
      out << reset_line << status.str();
      cout << out.str() << flush;
    }
  }
}

} // namespace
