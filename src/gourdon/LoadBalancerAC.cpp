///
/// @file  LoadBalancerAC.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
///        Load balancing is described in more detail at:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.md
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerAC.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace {

constexpr int64_t numbers_per_byte = primecount::SegmentedPiTable::numbers_per_byte();
constexpr int64_t l2_segment_size = L2_CACHE_SIZE * numbers_per_byte;

// Minimum segment size = 512 bytes.
// This size performs well on my AMD EPYC 2 near 1e16.
constexpr int64_t min_segment_size = (1 << 9) * numbers_per_byte;

} // namespace

namespace primecount {

LoadBalancerAC::LoadBalancerAC(int64_t sqrtx,
                               int64_t y,
                               int threads,
                               bool is_print) :
  sqrtx_(sqrtx),
  y_(y),
  threads_(threads),
  is_print_(is_print)
{
  lock_.init(threads);

  // The default segment size is x^(1/4). This
  // is tiny, will fit into the CPU's cache.
  int64_t x14 = isqrt(sqrtx);
  segment_size_ = x14;

  // When a single thread is used (and printing is
  // disabled) we can use a segment size larger
  // than x^(1/4) because load balancing is only
  // useful for multi-threading.
  if (threads == 1 && !is_print)
    segment_size_ = std::max(x14, l2_segment_size);

  segment_size_ = std::max(min_segment_size, segment_size_);
  segment_size_ = SegmentedPiTable::get_segment_size(segment_size_);

  // Most special leaves are below y (~ x^(1/3) * log(x)).
  // We make sure this interval is evenly distributed
  // amongst all threads by using a small segment size.
  // Above y we use a larger segment size but still ensure
  // that it fits into the CPU's cache.
  max_segment_size_ = std::max(l2_segment_size, segment_size_);

  if (is_print_)
    print_status(get_time());
}

bool LoadBalancerAC::get_work(ThreadDataAC& thread)
{
  double time = get_time();
  thread.secs = time - thread.secs;

  LockGuard lockGuard(lock_);

  if (low_ >= sqrtx_)
    return false;
  if (low_ == 0)
    start_time_ = time;

  double total_secs = time - start_time_;
  double increase_threshold = std::max(0.01, total_secs / 1000);
  int64_t remaining_dist = sqrtx_ - low_;

  // Most special leaves are below y (~ x^(1/3) * log(x)).
  // We make sure this interval is evenly distributed
  // amongst all threads by using a small segment size.
  // Above y we increase the segment size (or the number of
  // segments) by 2x if the thread runtime is close to 0.
  if (low_ > y_ &&
      thread.secs < increase_threshold &&
      thread.segment_size == segment_size_ &&
      segments_ * segment_size_ * (threads_ * 8) < remaining_dist)
  {
    int64_t increase_factor = 2;

    if (segment_size_ >= max_segment_size_)
      segments_ *= increase_factor;
    else
    {
      segment_size_ = segment_size_ * increase_factor;
      segment_size_ = std::min(segment_size_, max_segment_size_);
      segment_size_ = SegmentedPiTable::get_segment_size(segment_size_);
    }
  }

  if (is_print_)
    print_status(time);

  // Update current time because the thread may
  // have waited to acquire the lock.
  time = get_time();

  thread.low = low_;
  thread.segments = segments_;
  thread.segment_size = segment_size_;
  thread.secs = time;
  low_ = std::min(low_ + segments_ * segment_size_, sqrtx_);
  segment_nr_++;

  return thread.low < sqrtx_;
}

void LoadBalancerAC::print_status(double time)
{
  double threshold = 0.1;

  if (time - print_time_ >= threshold)
  {
    print_time_ = time;

    int64_t remaining_dist = sqrtx_ - low_;
    int64_t thread_dist = segments_ * segment_size_;
    int64_t total_segments = ceil_div(remaining_dist, thread_dist);
    total_segments += segment_nr_;

    std::ostringstream status;
    // Clear line because total_segments may become smaller
    status << "\r                                    "
            << "\rSegments: " << segment_nr_ << '/' << total_segments;
    std::cout << status.str() << std::flush;
  }
}

} // namespace
