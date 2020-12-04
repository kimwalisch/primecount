///
/// @file  StatusAC.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef STATUSAC_HPP
#define STATUSAC_HPP

#include <int128_t.hpp>

namespace primecount {

class StatusAC
{
public:
  StatusAC(maxint_t x);
  void print(int64_t b, int64_t max_b);
  double getPercent(int64_t b, int64_t max_b);
  void setPercent(double percent) { percent_ = percent; }
  void next();
private:
  bool isPrint(double time);
  void print(double percent);
  double epsilon_;
  double percent_ = -1;
  double percent_total_ = 0;
  double percent_segment_ = 80;
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double is_print_ = 0.1;
  int precision_;
};

} // namespace

#endif
