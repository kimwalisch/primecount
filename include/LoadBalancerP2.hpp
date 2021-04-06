///
/// @file  LoadBalancerP2.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERP2_HPP
#define LOADBALANCERP2_HPP

#include <stdint.h>

namespace primecount {

class LoadBalancerP2
{
public:
  LoadBalancerP2(int64_t low, int64_t z, int threads);
  int64_t get_thread_dist(int64_t low);
  int get_threads() const;

private:
  int64_t z_;
  int64_t min_dist_;
  int64_t thread_dist_;
  double time_;
  int threads_;
};

} // namespace

#endif
