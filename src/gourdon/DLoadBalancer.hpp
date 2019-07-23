///
/// @file  DLoadBalancer.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef DLOADBALANCER_HPP
#define DLOADBALANCER_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <LoadBalancer.hpp>

#include <stdint.h>

namespace primecount {

class DLoadBalancer
{
public:
  DLoadBalancer(maxint_t x, int64_t y, int64_t z);
  bool get_work(int64_t* low, int64_t* segments, int64_t* segment_size, maxint_t S2, Runtime& runtime);
  maxint_t get_sum() const;

private:
  void init_size();
  void update(int64_t* low, int64_t* segments, Runtime& runtime);
  void update_segments(Runtime& runtime);

  maxint_t sum_;
  int64_t low_;
  int64_t max_low_;
  int64_t xz_;
  int64_t segments_;
  int64_t segment_size_;
  int64_t max_size_;
  int64_t smallest_hard_leaf_;
  double time_;
};

} // namespace

#endif
