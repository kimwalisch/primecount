///
/// @file  StatusS2.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef STATUSS2_HPP
#define STATUSS2_HPP

#include <int128_t.hpp>

namespace primecount {

class StatusS2
{
public:
  StatusS2(maxint_t x, int64_t y, bool is_print);
  void print(double percent) const;
  void print(int64_t b, int64_t max_b);
  double getPercent(int64_t low, int64_t limit) const;
  double getStatus(int64_t low, int64_t limit);

private:
  int64_t y_log_y_ = 0;
  double epsilon_ = 0;
  double percent_ = -1;
  double x_tune_ = 0;
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double threshold_ = 0.1;
  int precision_ = 0;
};

} // namespace

#endif
