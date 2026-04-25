///
/// @file  LoadBalancerAC.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
///        Load balancing is described in more detail at:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERAC_HPP
#define LOADBALANCERAC_HPP

#include <primecount-config.hpp>
#include <macros.hpp>

#include <stdint.h>
#include <atomic>

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
  void store_packed(int64_t segment_size, int64_t segments);
  void print_AC_status(int64_t low, double time);

  int64_t sqrtx_ = 0;
  int64_t y_ = 0;
  int64_t min_segment_size_ = 0;
  int64_t max_segment_size_ = 0;
  double start_time_ = 0;
  int threads_ = 0;
  bool is_print_ = false;

  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<int64_t> low_{0};
  std::atomic<uint64_t> segment_data_{0};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  std::atomic<double> next_print_time_{0};
  std::atomic<bool> print_lock_{false};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

} // namespace

#endif
