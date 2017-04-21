///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent, requires
///        --status[=N] command-line flag.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace primecount {

S2Status::S2Status(maxint_t x) :
  old_percent_(0),
  old_time_(0),
  is_print_(1.0 / 20)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

bool S2Status::is_print(double time) const
{
  if (old_time_ == 0)
    return true;

  return (time - old_time_) >= is_print_;
}

double S2Status::skewed_percent(maxint_t n, maxint_t limit) const
{
  double exp = 0.96;
  double percent = get_percent((double) n, (double) limit);
  double base = exp + percent / (101 / (1 - exp));
  double low = pow(base, 100.0);
  percent = 100 - 100 * (pow(base, percent) - low) / (1 - low);

  return max(old_percent_, percent);
}

void S2Status::print(maxint_t n, maxint_t limit)
{
  double time = get_wtime();

  if (is_print(time))
  {
    #pragma omp critical (time)
    old_time_ = time;

    double percent = skewed_percent(n, limit);

    if ((percent > old_percent_) &&
        (percent - old_percent_) >= epsilon_)
    {
      #pragma omp critical (percent)
      old_percent_ = percent;

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

void S2Status::print(maxint_t n, maxint_t limit, double rsd)
{
  double time = get_wtime();

  if (is_print(time))
  {
    double percent = skewed_percent(n, limit);
    int load_balance = (int) in_between(0, 100 - rsd + 0.5, 100);

    ostringstream status;
    ostringstream out;

    status << "Status: " << fixed << setprecision(precision_) << percent << "%, ";
    status << "Load balance: " << load_balance << "%";
    size_t spaces = status.str().length() + 2;
    string reset_line = "\r" + string(spaces,' ') + "\r";
    out << reset_line << status.str();
    cout << out.str() << flush;

    old_percent_ = percent;
    old_time_ = time;
  }
}

} //namespace
