///
/// @file  Status.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef STATUS_HPP
#define STATUS_HPP

#include <int128_t.hpp>
#include <noinline.hpp>

namespace primecount {

class Status
{
public:
  Status(maxint_t x);
  NOINLINE void print(int64_t b, int64_t max_b);
  NOINLINE void print(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx);
  static double getPercent(int64_t low, int64_t limit, maxint_t sum, maxint_t sum_approx);
private:
  bool isPrint(double time);
  void print(double percent);
  double epsilon_;
  double percent_ = -1;
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double is_print_ = 0.1;
  int precision_;
};

} // namespace

#endif
