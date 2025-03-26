///
/// @file  LoadBalancerS2.hpp
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERS2_HPP
#define LOADBALANCERS2_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <OmpLock.hpp>
#include <StatusS2.hpp>

#include <stdint.h>

namespace primecount {

struct ThreadData
{
  int64_t low = 0;
  int64_t segments = 0;
  int64_t segment_size = 0;
  maxint_t sum = 0;
  double init_secs = 0;
  double secs = 0;

  void start_time()
  {
    secs = get_time();
  }

  void init_finished()
  {
    // Ensure start_time() has been called
    ASSERT(secs > 0);
    init_secs = get_time() - secs;
    ASSERT(init_secs >= 0);
  }

  void stop_time()
  {
    // Ensure start_time() has been called
    ASSERT(secs > 0);
    secs = get_time() - secs;
    ASSERT(secs >= 0);
  }
};

class LoadBalancerS2
{
public:
  LoadBalancerS2(maxint_t x, int64_t sieve_limit, maxint_t sum_approx, int threads, bool is_print);
  bool get_work(ThreadData& thread);
  maxint_t get_sum() const;

private:
  void update_load_balancing(const ThreadData& thread);
  void update_number_of_segments(const ThreadData& thread);
  double remaining_secs() const;

  int64_t low_ = 0;
  int64_t max_low_ = 0;
  int64_t sieve_limit_ = 0;
  int64_t sqrt_limit_ = 0;
  int64_t segments_ = 0;
  int64_t segment_size_ = 0;
  maxint_t sum_ = 0;
  maxint_t sum_approx_ = 0;
  double time_ = 0;
  int threads_ = 0;
  bool is_print_ = false;
  StatusS2 status_;
  OmpLock lock_;
};

} // namespace

#endif
