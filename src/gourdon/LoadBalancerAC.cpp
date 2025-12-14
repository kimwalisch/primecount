///
/// @file  LoadBalancerAC.cpp
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

#include "LoadBalancerAC.hpp"
#include "SegmentedPiTable.hpp"

#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>

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
  int64_t x14 = isqrt(sqrtx);

  // Minimum segment size = 512 bytes.
  // This size performs well near 1e16 on my AMD EPYC 2.
  int64_t min_segment_size = (1 << 9) * SegmentedPiTable::numbers_per_byte();

  // The maximum segment size matches the CPU's L1 cache
  // size (unless x^(1/4) > L1 cache size). This way
  // we ensure that most memory accesses will be cache
  // hits and we get good performance.
  int64_t L1_segment_size = L1_CACHE_SIZE * SegmentedPiTable::numbers_per_byte();

  if (threads == 1 && !is_print)
  {
    // When using a single thread (and printing is disabled)
    // we can use a segment size larger than x^(1/4)
    // because load balancing is only needed for multi-threading.
    segment_size_ = std::max(x14, L1_segment_size);
    segments_ = ceil_div(sqrtx, segment_size_);
  }
  else
  {
    // When using multi-threading we use a tiny segment size
    // of x^(1/4). This segment fits into the CPU's cache
    // and ensures good load balancing i.e. the work is evenly
    // distributed amongst all CPU cores.
    segment_size_ = x14;
    segments_ = 1;
  }

  segment_size_ = std::max(min_segment_size, segment_size_);
  segment_size_ = SegmentedPiTable::align_segment_size(segment_size_);
  max_segment_size_ = std::max(L1_segment_size, segment_size_);
  max_segment_size_ = SegmentedPiTable::align_segment_size(max_segment_size_);

  if (is_print_)
  {
    double time = get_time();
    std::string status = get_status(time);
    if (!status.empty())
      std::cout << status << std::flush;
  }
}

bool LoadBalancerAC::get_work(ThreadDataAC& thread)
{
  double time = get_time();
  thread.secs = time - thread.secs;
  std::string status;
  bool has_work;

  {
    LockGuard lockGuard(lock_);

    if (low_ >= sqrtx_)
      return false;
    if (low_ == 0)
      start_time_ = time;

    int64_t remaining_dist = sqrtx_ - low_;
    double total_secs = time - start_time_;
    double increase_threshold = std::max(0.01, total_secs / 1000);

    // Near the end of the computation we use a smaller
    // increase_threshold <= 1 second in order to make sure
    // all threads finish nearly at same time.
    if (segment_size_ == max_segment_size_)
      increase_threshold = std::min(increase_threshold, 1.0);

    // Most special leaves are below y (~ x^(1/3) * log(x)).
    // We make sure this interval is evenly distributed
    // amongst all threads by using a small segment size.
    // Above y we increase the segment size (or the number of
    // segments) by 2x if the thread runtime is close to 0.
    if (low_ > y_ &&
        thread.secs < increase_threshold &&
        thread.segments == segments_ &&
        thread.segment_size == segment_size_ &&
        segment_size_ * segments_ * (threads_ * 8) < remaining_dist)
    {
      int64_t increase_factor = 2;

      if (segment_size_ >= max_segment_size_)
        segments_ *= increase_factor;
      else
      {
        segment_size_ = segment_size_ * increase_factor;
        segment_size_ = std::min(segment_size_, max_segment_size_);
        segment_size_ = SegmentedPiTable::align_segment_size(segment_size_);
      }
    }

    if (is_print_)
      status = get_status(time);

    thread.low = low_;
    thread.segments = segments_;
    thread.segment_size = segment_size_;
    has_work = thread.low < sqrtx_;
    low_ = std::min(low_ + segment_size_ * segments_, sqrtx_);
    segment_nr_++;
  }

  // Printing to the terminal incurs a system call
  // and may hence be slow. Therefore, we do it
  // after having released the mutex.
  if (!status.empty())
    std::cout << status << std::flush;

  return has_work;
}

std::string LoadBalancerAC::get_status(double time)
{
  double threshold = 0.1;

  if (time - print_time_ >= threshold)
  {
    print_time_ = time;

    int64_t remaining_dist = sqrtx_ - low_;
    int64_t thread_dist = segment_size_ * segments_;
    int64_t total_segments = ceil_div(remaining_dist, thread_dist);
    total_segments += segment_nr_;

    std::string label = "Segments: ";
    std::string total_segs = std::to_string(total_segments);

    // Count characters in e.g. "Status: 1234/1234"
    std::size_t status_size = label.size() + total_segs.size() * 2 + 1;
    max_status_size_ = std::max(status_size, max_status_size_);
    std::string clear_line = '\r' + std::string(max_status_size_, ' ') + '\r';

    // The first part of the status string clears
    // the previous status line. This is necessary
    // because near the end of the computation the
    // status string becomes shorter.
    std::string status = clear_line + label;
    status += std::to_string(segment_nr_) + '/';
    status += total_segs;

    return status;
  }

  return std::string();
}

} // namespace
