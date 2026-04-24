///
/// @file  LoadBalancerAC.cpp
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

#include "LoadBalancerAC.hpp"
#include "SegmentedPiTable.hpp"

#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <print.hpp>
#include <TryLockGuard.hpp>

#include <stdint.h>
#include <atomic>
#include <string>

namespace primecount {

LoadBalancerAC::LoadBalancerAC(int64_t sqrtx,
                               int64_t y,
                               int threads,
                               bool is_print) :
  sqrtx_(sqrtx),
  y_(y),
  start_time_(get_time()),
  threads_(threads),
  is_print_(is_print)
{
  // When using multi-threading we use a tiny segment
  // size of x^(1/4). This segment fits into the CPU's
  // cache and ensures good load balancing i.e. the
  // work is evenly distributed amongst all CPU cores.
  int64_t x14 = isqrt(sqrtx);
  int64_t segment_size = x14;
  int64_t segments = 1;

  // Minimum segment size = 512 bytes.
  // This size performs well near 1e16 on my AMD EPYC 2.
  int64_t min_segment_size = (1 << 9) * SegmentedPiTable::numbers_per_byte();

  // The maximum segment size matches the CPU's L1 cache
  // size (unless x^(1/4) > L1 cache size). This way
  // we ensure that most memory accesses will be cache
  // hits and we get good performance.
  int64_t L1_segment_size = L1_CACHE_SIZE * SegmentedPiTable::numbers_per_byte();

  if (threads == 1)
  {
    // When using a single thread we can use a segment
    // size larger than x^(1/4) because load balancing
    // is only needed for multi-threading.
    segment_size = max(x14, L1_segment_size);

    if (!is_print)
      segments = ceil_div(sqrtx, segment_size);
  }

  // Limit to 2^31-1 to allow bit packing into uint64_t
  segments = min(segments, INT32_MAX);
  segment_size = min(segment_size, INT32_MAX);
  segment_size = max(min_segment_size, segment_size);
  segment_size = SegmentedPiTable::align_segment_size(segment_size);
  max_segment_size_ = max(L1_segment_size, segment_size);
  max_segment_size_ = SegmentedPiTable::align_segment_size(max_segment_size_);

  store_packed(segment_size, segments);

  if (is_print_)
  {
    ThreadDataAC thread;
    thread.low = 0;
    thread.segment_size = segment_size;
    thread.segments = segments;
    thread.secs = 0;
    print_AC_status(thread, start_time_);
  }
}

/// Pack segment_size & segments into a uint64_t,
/// needed for lockfree atomic data access.
///
void LoadBalancerAC::store_packed(int64_t segment_size,
                                  int64_t segments)
{
  ASSERT(segments <= INT32_MAX);
  ASSERT(segment_size <= INT32_MAX);
  uint64_t packed = segment_size | (segments << 32);
  segment_data_.store(packed, std::memory_order_relaxed);
}

/// Assign new [low, high[ workload to thread.
/// Multiple threads may call get_work() simultaneously, since
/// this function is not protected by a mutex, it must not
/// modify any shared member variables, except the atomic low_
/// and segment_data_ variables.
///
bool LoadBalancerAC::get_work(ThreadDataAC& thread)
{
  double time = get_time();
  thread.secs = time - thread.secs;
  int64_t low = low_.load(std::memory_order_relaxed);

  if (low >= sqrtx_)
    return false;

  int64_t remaining_dist = sqrtx_ - low;
  double total_secs = time - start_time_;
  double increase_threshold = max(0.01, total_secs / 1000);
  uint64_t segment_data = segment_data_.load(std::memory_order_relaxed);
  int64_t segment_size = segment_data & 0xffffffffu;
  int64_t segments = segment_data >> 32;

  // Near the end of the computation we use a smaller
  // increase_threshold <= 1 second in order to make sure
  // all threads finish nearly at same time.
  if (segment_size == max_segment_size_)
    increase_threshold = min(increase_threshold, 1.0);

  // Most special leaves are below y (~ x^(1/3) * log(x)).
  // We make sure this interval is evenly distributed
  // amongst all threads by using a small segment size.
  // Above y we increase the segment size (or the number of
  // segments) by 2x if the thread runtime is close to 0.
  if (low > y_ &&
      thread.secs < increase_threshold &&
      thread.segments >= segments &&
      thread.segment_size >= segment_size &&
      segment_size * segments * (threads_ * 8) < remaining_dist)
  {
    int64_t increase_factor = 2;

    if (segment_size >= max_segment_size_)
    {
      segments *= increase_factor;
      segments = min(segments, INT32_MAX);
    }
    else
    {
      segment_size = segment_size * increase_factor;
      segment_size = min(segment_size, max_segment_size_);
      segment_size = SegmentedPiTable::align_segment_size(segment_size);
    }

    store_packed(segment_size, segments);
  }

  int64_t thread_dist = segment_size * segments;
  low = low_.fetch_add(thread_dist, std::memory_order_relaxed);

  thread.low = low;
  thread.segments = segments;
  thread.segment_size = segment_size;

  // The lockfree critical section above should complete
  // as fast as possible. Hence, printing should be done
  // afterwards since it may incur a system call.
  if (is_print_)
    print_AC_status(thread, time);

  return thread.low < sqrtx_;
}

void LoadBalancerAC::print_AC_status(const ThreadDataAC& thread,
                                     double time)
{
  int64_t segment_nr = segment_nr_.fetch_add(1, std::memory_order_relaxed);

#if __cplusplus >= 201703L
  // Prevent lock contention on many-core systems,
  // print status only every 0.1 seconds.
  if (std::atomic<double>::is_always_lock_free &&
      time <= next_print_time_.load(std::memory_order_relaxed))
    return;
#endif

  // For printing the status it is OK to use a non-blocking
  // userspace lock because printing the status is a non
  // essential operation and hence even if the OS preempts
  // the thread holding the lock it won't cause any deadlocks
  // or performance issues, it will only delay the status
  // output.
  TryLockGuard guard(print_lock_);

  if (guard.owns_lock())
  {
  #if __cplusplus >= 201703L
    // It is theoretically possible that multiple threads
    // enter this critical sections within 0.1 seconds.
    // This additional condition prevents it.
    if (std::atomic<double>::is_always_lock_free &&
        time <= next_print_time_.load(std::memory_order_relaxed))
      return;
  #endif

    // The next thread can print again in 0.1 seconds
    next_print_time_.store(time + 0.1, std::memory_order_relaxed);

    segment_nr += 1;
    int64_t remaining_dist = sqrtx_ - thread.low;
    int64_t thread_dist = thread.segment_size * thread.segments;
    int64_t total_segments = ceil_div(remaining_dist, thread_dist);
    total_segments += segment_nr;

    std::string status = "Segments: ";
    status += std::to_string(segment_nr) + '/';
    status += std::to_string(total_segments);
    print_status(status);
  }
}

} // namespace
