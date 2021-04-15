///
/// @file  LoadBalancerAC.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
///        Load balancing is described in more detail at:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.md
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERAC_HPP
#define LOADBALANCERAC_HPP

#include <OmpLock.hpp>
#include <stdint.h>

namespace primecount {

class LoadBalancerAC
{
public:
  LoadBalancerAC(int64_t sqrtx, int64_t y, int threads, bool is_print);
  bool get_work(int64_t& low, int64_t& high);

private:
  void validate_segment_sizes();
  void compute_total_segments();
  void print_status();

  int64_t low_ = 0;
  int64_t sqrtx_ = 0;
  int64_t x14_ = 0;
  int64_t y_ = 0;
  int64_t segment_size_ = 0;
  int64_t large_segment_size_ = 0;
  int64_t segment_nr_ = 0;
  int64_t total_segments_ = 0;
  double time_ = 0;
  int threads_ = 0;
  bool is_print_ = false;
  OmpLock lock_;
};

} // namespace

#endif
