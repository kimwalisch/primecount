///
/// @file  MpiLoadBalancer.cpp
/// @brief The MpiLoadBalancer evenly distributes the computation
///        of the hard special leaves onto cluster nodes.
///
///        Simply parallelizing the computation of the special
///        leaves in the Lagarias-Miller-Odlyzko algorithm by
///        subdividing the sieve interval by the number of threads
///        into equally sized subintervals does not scale because
///        the distribution of the special leaves is highly skewed
///        and most special leaves are in the first few segments
///        whereas later on there are very few special leaves.
///
///        This MpiLoadBalancer gradually increases the number of
///        segments to sieve as long the expected runtime of the
///        sieve distance is smaller than the expected finish time
///        of the algorithm. Near the end the MpiLoadBalancer will
///        gradually decrease the number of segments to sieve in
///        order to prevent that 1 thread will run much longer
///        than all the other threads.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <MpiLoadBalancer.hpp>
#include <MpiMsg.hpp>
#include <Status.hpp>
#include <Sieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;

namespace primecount {

MpiLoadBalancer::MpiLoadBalancer(maxint_t x,
                                 int64_t sieve_limit,
                                 maxint_t sum_approx) :
  low_(0),
  max_low_(0),
  sieve_limit_(sieve_limit),
  segments_(1),
  sum_(0),
  sum_approx_(sum_approx),
  time_(get_time()),
  status_(x)
{
  // start with a tiny segment_size as most
  // special leaves are in the first few segments
  // and we need to ensure that all threads are
  // assigned an equal amount of work
  int64_t sqrt_limit = isqrt(sieve_limit);
  int64_t log = ilog(sqrt_limit);
  log = max(log, 1);
  segment_size_ = sqrt_limit / log;

  int64_t min_size = 1 << 9;
  segment_size_ = max(segment_size_, min_size);
  segment_size_ = Sieve::get_segment_size(segment_size_);

  // try to use a segment size that fits exactly
  // into the CPUs L1 data cache
  int64_t l1_dcache_size = 1 << 15;
  max_size_ = l1_dcache_size * 30;
  max_size_ = max(max_size_, sqrt_limit);
  max_size_ = Sieve::get_segment_size(max_size_);
}

void MpiLoadBalancer::get_work(MpiMsg* msg)
{
  sum_ += msg->sum();

  if (msg->low() > max_low_)
  {
    max_low_ = msg->low();
    segments_ = msg->segments();
    segment_size_ = msg->segment_size();

    // We only start increasing the segment_size and segments
    // per thread once the first special leaves have been
    // found. Near the start there is a very large number of
    // leaves and we don't want a single thread to compute
    // them all by himself (which would cause scaling issues).
    if (sum_ != 0)
    {
      if (segment_size_ < max_size_)
        segment_size_ = min(segment_size_ * 2, max_size_);
      else
      {
        ThreadSettings thread;
        thread.init_secs = msg->init_seconds();
        thread.secs = msg->seconds();
        update_segments(thread);
      }
    }
  }

  // udpate msg with new work todo
  msg->update(low_, segments_, segment_size_);
  low_ += segments_ * segment_size_;
  low_ = min(low_, sieve_limit_);
}

/// Remaining seconds till finished
double MpiLoadBalancer::remaining_secs() const
{
  double percent = status_.getPercent(low_, sieve_limit_, sum_, sum_approx_);
  percent = in_between(10, percent, 100);
  double total_secs = get_time() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

/// Increase or decrease the number of segments based on the
/// remaining runtime. Near the end it is important that
/// threads run only for a short amount of time in order to
/// ensure all threads finish nearly at the same time.
///
void MpiLoadBalancer::update_segments(ThreadSettings& thread)
{
  double rem = remaining_secs();
  double threshold = rem / 3;
  double min_secs = 0.01;

  // Each thread should run at least 10x
  // longer than its initialization time
  threshold = max(threshold, thread.init_secs * 10);
  threshold = max(threshold, min_secs);

  // divider must not be 0
  double divider = max(thread.secs, min_secs / 10);
  double factor = threshold / divider;

  // Reduce the thread runtime if it is much
  // larger than its initialization time
  if (thread.secs > min_secs &&
      thread.secs > thread.init_secs * 1000)
  {
    double old = factor;
    factor = (thread.init_secs * 1000) / thread.secs;
    factor = min(factor, old);
  }

  factor = in_between(0.5, factor, 2.0);
  double new_segments = round(segments_ * factor);
  segments_ = (int64_t) new_segments;
  segments_ = max(segments_, 1);
}

} // namespace
