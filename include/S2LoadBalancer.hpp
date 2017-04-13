///
/// @file  S2LoadBalancer.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2LOADBALANCER_HPP
#define S2LOADBALANCER_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <stdint.h>

namespace primecount {

class S2LoadBalancer
{
public:
  S2LoadBalancer(maxint_t x, int64_t y, int64_t z, int64_t threads);
  S2LoadBalancer(maxint_t x, int64_t y, int64_t z, int64_t threads, double rsd);
  int64_t get_min_segment_size() const;
  double get_rsd() const;
  void update(int64_t* segment_size,
              int64_t* segments_per_thread,
              int64_t low,
              int64_t threads,
              thread_timings_t& timings);
private:
  void init(maxint_t x, int64_t y, int64_t threads);
  void set_min_size(int64_t z);
  void update(int64_t* segments_per_thread, double seconds, double pivot);
  void update_min_size(double divisor);
  double get_avg_seconds() const;
  double get_pivot(double seconds) const;
  bool is_increase(double seconds, double pivot) const;
  bool is_decrease(double seconds, double pivot) const;
  double x_;
  double y_;
  double z_;
  double rsd_;
  double count_;
  double total_seconds_;
  double min_seconds_;
  double decrease_dividend_;
  int64_t min_size_;
  int64_t sqrtz_;
  int64_t smallest_hard_leaf_;
};

} // namespace

#endif
