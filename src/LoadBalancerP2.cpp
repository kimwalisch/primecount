///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerP2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <print.hpp>
#include <TryLockGuard.hpp>

#include <stdint.h>
#include <atomic>
#include <cmath>
#include <string>

namespace primecount {

/// We need to sieve [sqrt(x), sieve_limit[
LoadBalancerP2::LoadBalancerP2(maxint_t x,
                               int64_t sieve_limit,
                               int threads,
                               bool is_print) :
  sieve_limit_(sieve_limit),
  precision_(get_status_precision(x)),
  is_print_(is_print)
{
  int64_t low = isqrt(x);
  low = min(low, sieve_limit_);
  low_.store(low, std::memory_order_relaxed);
  int64_t dist = sieve_limit_ - low;

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  min_thread_dist_ = 1 << 23;
  double max_threads = std::pow(sieve_limit_, 1 / 3.7);
  threads = min(threads, (int) max_threads);
  threads_ = ideal_num_threads(dist, threads, min_thread_dist_);

  // Using more chunks per thread improves load
  // balancing but also adds some overhead.
  int64_t chunks_per_thread = 8;
  thread_dist_ = dist / (threads_ * chunks_per_thread);
  thread_dist_ = max(min_thread_dist_, thread_dist_);
}

int LoadBalancerP2::get_threads() const
{
  return threads_;
}

/// Assign new [low, high[ workload to thread.
/// Multiple threads may call get_work() simultaneously, since
/// this function is not protected by a mutex, it must not
/// modify any shared member variables, except the atomic low_.
///
bool LoadBalancerP2::get_work(int64_t& low, int64_t& high)
{
  low = low_.load(std::memory_order_relaxed);

  if (low >= sieve_limit_)
    return false;

  int64_t dist = sieve_limit_ - low;
  int64_t thread_dist = thread_dist_;

  if (threads_ == 1)
  {
    // When using a single thread (and printing is disabled) we
    // can set thread_dist to the entire sieving distance since
    // load balancing is only useful for multi-threading.
    if (!is_print_)
      thread_dist = dist;
  }
  else
  {
    // Ensure that the thread initialization i.e. the calculation of
    // PrimePi(low) uses less time than the actual computation.
    // Computing PrimePi(low) uses O(low^(2/3) / log(low)^2) time
    // whereas sieving a distance of n = low^(2/3) uses O(n log log n)
    // time. Hence, the sieving time is larger than the initialization
    // time. Using these settings for low = 2e14 the sieving time is
    // 5x larger than the initialization time on my AMD EPYC 7642 CPU.
    double low13 = std::cbrt(low);
    int64_t low23 = (int64_t) (low13 * low13);
    int64_t min_thread_dist = max(min_thread_dist_, low23);
    thread_dist = max(min_thread_dist, thread_dist);

    // Reduce thread distance near the end to keep all
    // threads busy until the computation finishes.
    int64_t max_thread_dist = dist / threads_;
    if (thread_dist > max_thread_dist)
      thread_dist = max(min_thread_dist, max_thread_dist);
  }

  // The earlier loads are used for heuristic chunk
  // sizing, it is OK if they are slightly outdated. This
  // fetch_add() reserves unique work for this thread.
  low = low_.fetch_add(thread_dist, std::memory_order_relaxed);
  high = low + thread_dist;
  high = min(high, sieve_limit_);
  bool has_work = low < sieve_limit_;

  // The lockfree critical section above should complete
  // as fast as possible. Hence, printing should be done
  // afterwards since it may incur a system call.
  if (is_print_ && has_work)
    print_P2_status(low);

  return has_work;
}

void LoadBalancerP2::print_P2_status(int64_t low)
{
  double time = get_time();

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
    double percent = get_percent(low, sieve_limit_);

    if (percent > percent_)
    {
      percent_ = percent;
      std::string status = "Status: ";
      status += to_string(percent, precision_);
      status += '%';
      print_status(status);
    }
  }
}

} // namespace
