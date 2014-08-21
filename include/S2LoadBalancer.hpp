///
/// @file  S2LoadBalancer.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2LOADBALANCER_HPP
#define S2LOADBALANCER_HPP

#include <aligned_vector.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

class S2LoadBalancer
{
public:
  S2LoadBalancer(maxint_t x, int64_t sqrtz, int64_t threads);
  int64_t get_min_segment_size() const;
  double get_rsd() const;
  void update(int64_t low,
              int64_t threads,
              int64_t* segment_size,
              int64_t* segments_per_thread,
              aligned_vector<double>& timings);
private:
  void set_min_size(int64_t z);
  void update_avg_seconds(double seconds);
  void update_min_size(double divisor);
  double get_decrease_threshold(double seconds, int64_t threads) const;
  bool increase_size(double seconds, double decrease) const;
  bool decrease_size(double seconds, double decrease) const;
  double x_;
  double z_;
  double rsd_;
  double avg_seconds_;
  double min_seconds_;
  int64_t min_size_;
  int64_t max_size_;
  int64_t count_;
};

} // namespace

#endif
