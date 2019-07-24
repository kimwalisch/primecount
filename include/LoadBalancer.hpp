///
/// @file  LoadBalancer.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCER_HPP
#define LOADBALANCER_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <S2Status.hpp>

#include <stdint.h>

namespace primecount {

struct Runtime
{
  Runtime() { reset(); }
  void reset() { init = 0; secs = 0; }
  void start() { reset(); secs = get_time(); }
  void stop() { secs = get_time() - secs; }
  void init_start() { init = get_time(); }
  void init_stop() { init = get_time() - init; }
  double init;
  double secs;
};

class LoadBalancer
{
public:
  LoadBalancer(maxint_t x, int64_t sieve_limit, maxint_t sum_approx);
  bool get_work(int64_t* low, int64_t* segments, int64_t* segment_size, maxint_t sum, Runtime& runtime);
  maxint_t get_sum() const;

private:
  void update(int64_t* low, int64_t* segments, Runtime& runtime);
  void update_segments(Runtime& runtime);
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
  S2Status status_;
};

} // namespace

#endif
