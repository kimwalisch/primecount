///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerP2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;

namespace primecount {

LoadBalancerP2::LoadBalancerP2(int64_t low,
                               int64_t z,
                               int threads) :
  z_(z),
  min_dist_(1 << 22),
  thread_dist_(min_dist_),
  time_(-1)
{
  int64_t dist = z_ - min(low, z_);
  threads_ = ideal_num_threads(threads, dist, thread_dist_);
}

int LoadBalancerP2::get_threads() const
{
  return threads_;
}

void LoadBalancerP2::set_thread_dist(int64_t thread_dist)
{
  thread_dist_ = thread_dist;
}

/// Calculate the thread sieving distance for the next
/// iteration. Since all threads must synchronize after
/// each iteration we want to gradually increase the
/// thread distance in order to ensure that all threads
/// run for approximately the same amount of time.
///
int64_t LoadBalancerP2::get_thread_dist(int64_t low)
{
  double start_time = time_;
  time_ = get_time();
  double seconds = time_ - start_time;

  if (start_time > 0)
  {
    if (seconds < 60)
      thread_dist_ *= 2;
    if (seconds > 60)
      thread_dist_ /= 2;
  }

  int64_t dist = z_ - min(low, z_);
  int64_t max_dist = ceil_div(dist, threads_);
  thread_dist_ = in_between(min_dist_, thread_dist_, max_dist);
  int64_t next_dist = threads_ * thread_dist_;
  int64_t last_dist = threads_ * min_dist_;

  // Keep all threads busy in the last iteration 
  if (low + next_dist + last_dist > dist)
    thread_dist_ = max(min_dist_, max_dist);

  return thread_dist_;
}

} // namespace
