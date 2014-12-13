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

#include <inttypes.hpp>

namespace primecount {

class S2Status
{
public:
  S2Status();
  void print(maxint_t n, maxint_t limit);
  void print(maxint_t n, maxint_t limit, double rsd);
private:
  int calculate_percent(maxint_t n, maxint_t limit) const;
  int old_;
  double time_;
};

} // namespace

#endif
