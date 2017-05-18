///
/// @file  S2Status.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2STATUS_HPP
#define S2STATUS_HPP

#include <int128_t.hpp>
#include <atomic>

namespace primecount {

class S2Status
{
public:
  S2Status(maxint_t x);
  void print(maxint_t n, maxint_t limit);
  void print(maxint_t n, maxint_t limit, double rsd);
  static double getPercent(int64_t low, int64_t limit, maxint_t S2, maxint_t S2_approx);
private:
  bool is_print(double time) const;
  void print(double percent) const;
  static double skewed_percent(maxint_t x, maxint_t y);
  double epsilon_;
  std::atomic<double> percent_;
  std::atomic<double> time_;
  double is_print_;
  int precision_;
};

} // namespace

#endif
