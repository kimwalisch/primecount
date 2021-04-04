///
/// @file  StatusAC.hpp
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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
  void print(int64_t low, int64_t limit, int64_t segment_size);
private:
  bool isPrint(double time);
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double is_print_ = 0.1;
};

} // namespace

#endif
