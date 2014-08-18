///
/// @file  S2Status.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primecount-internal.hpp>
#include <ptypes.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace primecount {

S2Status::S2Status(maxint_t s2_approx) :
  s2_approx_((double) s2_approx),
  percent_(0)
{ }

void S2Status::print(maxint_t s2_current, double rsd)
{
  double percent = get_percent((double) s2_current, s2_approx_);
  double base = 0.96 + percent / (101 / (1 - 0.96));
  double min = pow(base, 100.0);
  double max = pow(base, 0.0);
  percent = 100 - in_between(0, 100 * (pow(base, percent) - min) / (max - min), 100);
  percent_ = std::max(percent_, (int) percent);

  int load_balance = (int) in_between(0, 100 - rsd + 0.5, 100);

  ostringstream oss;
  oss << "\r" << string(40,' ') << "\r";
  oss << "Status: " << (int) percent_ << "%, ";
  oss << "Load balance: " << load_balance << "%";
  cout << oss.str() << flush;
}

} //namespace
