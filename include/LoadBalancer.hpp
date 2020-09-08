///
/// @file  LoadBalancer.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCER_HPP
#define LOADBALANCER_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <OmpLock.hpp>
#include <Status.hpp>
#include <noinline.hpp>

#include <stdint.h>

namespace primecount {

struct ThreadSettings
{
  int64_t low = 0;
  int64_t segments = 0;
  int64_t segment_size = 0;
  maxint_t sum = 0;
  double init_secs = 0;
  double secs = 0;

  void start_time() { secs = get_time(); }
  void stop_time() { secs = get_time() - secs; }
  void init_finished() { init_secs = get_time() - secs; }
};

class LoadBalancer
{
public:
  NOINLINE LoadBalancer(maxint_t x, int64_t sieve_limit, maxint_t sum_approx);
  NOINLINE bool get_work(ThreadSettings& thread);
  maxint_t get_sum() const;

private:
  void update(ThreadSettings& thread);
  void update_segments(ThreadSettings& thread);
  double remaining_secs() const;

  int64_t low_;
  int64_t max_low_;
  int64_t sieve_limit_;
  int64_t segments_;
  int64_t segment_size_;
  int64_t max_size_;
  maxint_t sum_;
  maxint_t sum_approx_;
  double time_;
  Status status_;
  OmpLock lock_;
};

} // namespace

#endif
