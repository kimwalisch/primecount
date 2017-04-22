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
  double skewed_percent(maxint_t n, maxint_t limit) const;
private:
  bool is_print(double time) const;
  double epsilon_;
  std::atomic<double> percent_;
  std::atomic<double> time_;
  double is_print_;
  int precision_;
};

} // namespace

#endif
