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
  StatusS2(maxint_t x, bool is_print);
  StatusS2(maxint_t x, double leaves_low, int64_t limit, double sieve_factor, bool is_print);
  double getPercent(int64_t low) const;
  void print_S2_hard(double percent);
  void print_S2_easy(int64_t b, int64_t max_b);

private:
  double work(double low) const;
  double leaves_low_ = 0;
  double leaves_scale_ = 1;
  double sieve_work_ = 0;
  double limit_ = 1;
  double total_work_ = 1;
  double epsilon_ = 0;
  double percent_ = -1;
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double threshold_ = 0.1;
  int precision_ = 0;
};

} // namespace

#endif
