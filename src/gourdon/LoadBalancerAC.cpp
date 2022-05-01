///
/// @file  LoadBalancerAC.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
///        Load balancing is described in more detail at:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.md
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerAC.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>

namespace {

// CPU L2 cache size per core (in bytes)
constexpr int64_t l2_cache_bytes = 512 << 10;
constexpr int64_t numbers_per_byte = primecount::SegmentedPiTable::numbers_per_byte();
constexpr int64_t l2_segment_size = l2_cache_bytes * numbers_per_byte;

// Minimum segment size = 1 KiB
constexpr int64_t min_segment_size = (1 << 10) * numbers_per_byte;

} // namespace

namespace primecount {

LoadBalancerAC::LoadBalancerAC(int64_t sqrtx,
                               int64_t y,
                               int threads,
                               bool is_print) :
  sqrtx_(sqrtx),
  x14_(isqrt(sqrtx)),
  y_(y),
  threads_(threads),
  is_print_(is_print)
{
  if (threads > 1)
    lock_.init();

  // When a single thread is used (and printing is
  // disabled) we can use a segment size larger
  // than x^(1/4) because load balancing is only
  // useful for multi-threading.
  if (threads == 1 && !is_print)
    segment_size_ = std::max(x14_, l2_segment_size);
  else
  {
    // The default segment size is x^(1/4). This
    // is tiny, will fit into the CPU's cache.
    segment_size_ = x14_;

    // Most special leaves are below y (~ x^(1/3) * log(x)).
    // We make sure this interval is evenly distributed
    // amongst all threads by using a small segment size.
    // Above y we use a larger segment size but still ensure
    // that it fits into the CPU's cache.
    if (y_ < sqrtx_)
    {
      int64_t max_segment_size = (sqrtx_ - y_) / (threads_ * 8);
      large_segment_size_ = segment_size_ * 16;
      large_segment_size_ = min3(large_segment_size_, l2_segment_size, max_segment_size);
      large_segment_size_ = std::max(segment_size_, large_segment_size_);
    }
  }

  validate_segment_sizes();
  compute_total_segments();
  print_status();
}

bool LoadBalancerAC::get_work(int64_t& low, int64_t& high)
{
  LockGuard lockGuard(lock_);

  if (low_ >= sqrtx_)
    return false;

  // Most special leaves are below y (~ x^(1/3) * log(x)).
  // We make sure this interval is evenly distributed
  // amongst all threads by using a small segment size.
  // Above y we use a larger segment size but still ensure
  // that it fits into the CPU's cache.
  if (low_ > y_)
    segment_size_ = large_segment_size_;

  low = low_;
  high = low + segment_size_;
  high = std::min(high, sqrtx_);
  low_ = high;
  segment_nr_++;
  print_status();

  return low < sqrtx_;
}

void LoadBalancerAC::validate_segment_sizes()
{
  segment_size_ = std::max(min_segment_size, segment_size_);
  large_segment_size_ = std::max(segment_size_, large_segment_size_);
  segment_size_ = SegmentedPiTable::get_segment_size(segment_size_);
  large_segment_size_ = SegmentedPiTable::get_segment_size(large_segment_size_);
}

void LoadBalancerAC::compute_total_segments()
{
  int64_t small_segments = ceil_div(y_, segment_size_);
  int64_t threshold = std::min(small_segments * segment_size_, sqrtx_);
  int64_t large_segments = ceil_div(sqrtx_ - threshold, large_segment_size_);
  total_segments_ = small_segments + large_segments;
}

void LoadBalancerAC::print_status()
{
  if (is_print_)
  {
    double time = get_time();
    double old = time_;
    double threshold = 0.1;

    if ((time - old) >= threshold)
    {
      time_ = time;
      std::cout << "\rSegments: " << segment_nr_ << "/" << total_segments_ << std::flush;
    }
  }
}

} // namespace
