///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERP2_HPP
#define LOADBALANCERP2_HPP

#include <noinline.hpp>
#include <stdint.h>

namespace primecount {

class LoadBalancerP2
{
public:
  NOINLINE LoadBalancerP2(int64_t z, int threads);
  NOINLINE int64_t get_thread_dist(int64_t low);
  int get_threads() const;

private:
  double time_ = -1;
  int64_t min_dist_ = 1 << 22;
  int64_t thread_dist_ = min_dist_;
  int threads_;
  int64_t z_;
};

} // namespace

#endif
