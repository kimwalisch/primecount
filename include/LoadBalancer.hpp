///
/// @file  LoadBalancer.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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
  void reset() { init = 0; seconds = 0; }
  void start() { reset(); seconds = get_wtime(); }
  void stop() { seconds = get_wtime() - seconds; }
  void init_start() { init = get_wtime(); }
  void init_stop() { init = get_wtime() - init; }
  double init;
  double seconds;
};

class LoadBalancer
{
public:
  LoadBalancer(maxint_t x, int64_t y, maxint_t s2_approx);
  void get_work(int64_t* low, int64_t* segments, int64_t* segment_size, maxint_t S2, Runtime& runtime);
  bool finished();
  maxint_t get_result();

private:
  void init_size(maxint_t x, int64_t y);
  bool is_increase(double percent, Runtime& runtime);

  S2Status status_;
  int64_t low_;
  int64_t limit_;
  int64_t segments_;
  int64_t segment_size_;
  int64_t max_size_;
  maxint_t s2_approx_;
  maxint_t S2_total_;
  double time_;
  bool finished_;
};

} // namespace

#endif
