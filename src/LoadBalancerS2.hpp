///
/// @file  LoadBalancerS2.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERS2_HPP
#define LOADBALANCERS2_HPP

#include <primecount-internal.hpp>
#include <primecount-config.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <StatusS2.hpp>

#include <stdint.h>
#include <atomic>

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
  LoadBalancerS2(maxint_t x, int64_t y, int64_t sieve_limit, int threads, bool is_print);
  bool get_work(ThreadData& thread);

private:
  double remaining_secs(int64_t low) const;
  void store_packed(uint64_t segment_size, uint64_t segments);
  void run_load_balancing(ThreadData& thread, int64_t segment_size);
  int64_t get_segments(const ThreadData& thread, int64_t low) const;
  void print_S2_status(int64_t high);

  int64_t sieve_limit_ = 0;
  int64_t sqrt_limit_ = 0;
  double start_time_ = 0;
  int threads_ = 0;
  bool is_print_ = false;
  StatusS2 status_;

  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<int64_t> low_{0};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  std::atomic<int64_t> max_low_{0};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
  std::atomic<uint64_t> segment_data_{0};
  MAYBE_UNUSED char pad4[MAX_CACHE_LINE_SIZE];
  std::atomic<double> next_print_time_{0};
  std::atomic<bool> print_lock_{false};
  std::atomic<bool> found_first_leaf_{false};
  MAYBE_UNUSED char pad5[MAX_CACHE_LINE_SIZE];
};

} // namespace

#endif
