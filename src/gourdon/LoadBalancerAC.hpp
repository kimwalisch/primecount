///
/// @file  LoadBalancerAC.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
///        Load balancing is described in more detail at:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERAC_HPP
#define LOADBALANCERAC_HPP

#include <OmpLock.hpp>

#include <stdint.h>
#include <string>

namespace primecount {

struct ThreadDataAC
{
  int64_t low = 0;
  int64_t segments = 0;
  int64_t segment_size = 0;
  double secs = 0;
};

class LoadBalancerAC
{
public:
  LoadBalancerAC(int64_t sqrtx, int64_t y, int threads, bool is_print);
  bool get_work(ThreadDataAC& thread);

private:
  std::string get_status(double current_time);
  int64_t low_ = 0;
  int64_t sqrtx_ = 0;
  int64_t y_ = 0;
  int64_t segments_ = 0;
  int64_t segment_size_ = 0;
  int64_t segment_nr_ = 0;
  int64_t max_segment_size_ = 0;
  double start_time_ = 0;
  double print_time_ = 0;
  int threads_ = 0;
  bool is_print_ = false;
  OmpLock lock_;
};

} // namespace

#endif
