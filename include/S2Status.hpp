///
/// @file  S2Status.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2STATUS_HPP
#define S2STATUS_HPP

#include <int128.hpp>

namespace primecount {

class S2Status
{
public:
  S2Status(maxint_t s2_approx);
  void print(maxint_t s2_current, double rsd);
private:
  double s2_approx_;
  int percent_;
};

} // namespace

#endif
