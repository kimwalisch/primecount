///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent.
///        Requires --status command-line flag.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primecount-internal.hpp>
#include <pmath.hpp>
#include <int128.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace primecount {

S2Status::S2Status() :
  old_(-1),
  time_(0)
{ }

int S2Status::skewed_percent(maxint_t n, maxint_t limit) const
{
  double exp = 0.96;
  double percent = get_percent((double) n, (double) limit);
  double base = exp + percent / (101 / (1 - exp));
  double min = pow(base, 100.0);
  percent = 100 - in_between(0, 100 * (pow(base, percent) - min) / (1 - min), 100);
  return max(old_, (int) percent);
}

void S2Status::print(maxint_t n, maxint_t limit)
{
  int percent = skewed_percent(n, limit);
  if (percent > old_)
  {
    #pragma omp critical (s2_status)
    old_ = percent;

    ostringstream oss;
    oss << "\r" << string(12,' ');
    oss << "\rStatus: " << percent << "%";
    cout << oss.str() << flush;
  }
}

void S2Status::print(maxint_t n, maxint_t limit, double rsd)
{
  double t2 = get_wtime();
  if (old_ >= 0 && (t2 - time_) < 0.01)
    return;

  time_ = t2;
  int percent = skewed_percent(n, limit);
  int load_balance = (int) in_between(0, 100 - rsd + 0.5, 100);
  old_ = percent;

  ostringstream oss;
  oss << "\r" << string(40,' ');
  oss << "\rStatus: " << percent << "%, ";
  oss << "Load balance: " << load_balance << "%";
  cout << oss.str() << flush;
}

} //namespace
