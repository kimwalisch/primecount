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

#include <int128_t.hpp>
#include <OmpLock.hpp>

#include <stdint.h>

namespace primecount {

class LoadBalancerP2
{
public:
  LoadBalancerP2(maxint_t x, int64_t sieve_limit, int threads, bool is_print);
  bool get_work(int64_t& low, int64_t& high);
  int get_threads() const;

private:
  void print_status();

  int64_t low_;
  int64_t sieve_limit_;
  int64_t thread_dist_;
  double precision_;
  int threads_;
  bool is_print_;
  OmpLock lock_;
};

} // namespace

#endif
