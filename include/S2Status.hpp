///
/// @file  S2Status.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2STATUS_HPP
#define S2STATUS_HPP

#include <int128_t.hpp>
#include <noinline.hpp>

#if defined(_OPENMP)
  #include <OmpLock.hpp>
#endif

namespace primecount {

class S2Status
{
public:
  S2Status(maxint_t x);
  NOINLINE void print(maxint_t n, maxint_t limit);
  void setPercent(double percent) { percent_ = percent; }
  static double getPercent(int64_t low, int64_t limit, maxint_t S2, maxint_t S2_approx);
private:
  bool is_print(double time);
  static double skewed_percent(maxint_t x, maxint_t y);
  double epsilon_;
  double percent_ = -1;
  double time_ = 0;
  double is_print_ = 1.0 / 20;
  int precision_;

#if defined(_OPENMP)
  OmpLock lock_;
#endif
};

} // namespace

#endif
