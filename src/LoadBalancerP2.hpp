///
/// @file  LoadBalancerP2.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERP2_HPP
#define LOADBALANCERP2_HPP

#include <primecount-config.hpp>
#include <int128_t.hpp>
#include <macros.hpp>

#include <stdint.h>
#include <atomic>

namespace primecount {

class LoadBalancerP2
{
public:
  LoadBalancerP2(maxint_t x, int64_t sieve_limit, int threads, bool is_print);
  bool get_work(int64_t& low, int64_t& high);
  int get_threads() const;

private:
  void print_P2_status(int64_t low);

  int64_t sieve_limit_ = 0;
  int64_t min_thread_dist_ = 0;
  int64_t thread_dist_ = 0;
  double percent_ = -1;
  int threads_ = 0;
  int precision_ = 0;
  bool is_print_ = false;

  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<int64_t> low_{0};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  std::atomic<double> next_print_time_{0};
  std::atomic<bool> print_lock_{false};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

} // namespace

#endif
