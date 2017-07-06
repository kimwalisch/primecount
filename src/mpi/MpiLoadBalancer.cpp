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
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <MpiLoadBalancer.hpp>
#include <MpiMsg.hpp>
#include <imath.hpp>
#include <Sieve.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;

namespace primecount {

MpiLoadBalancer::MpiLoadBalancer(maxint_t x,
                                 int64_t y,
                                 int64_t z,
                                 maxint_t s2_approx) :
  z_(z),
  limit_(z + 1),
  low_(0),
  max_low_(0),
  segments_(1),
  s2_hard_(0),
  s2_approx_(s2_approx),
  time_(get_wtime()),
  status_(x)
{
  double alpha = get_alpha(x, y);
  maxint_t x16 = iroot<6>(x);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * x16));

  // try to use a segment size that fits exactly
  // into the CPUs L1 data cache
  int64_t l1_dcache_size = 1 << 15;
  int64_t size = l1_dcache_size * 30;
  size = max(size, isqrt(z));
  segment_size_ = Sieve::get_segment_size(size);
}

void MpiLoadBalancer::get_work(MpiMsg* msg, maxint_t s2_hard)
{
  s2_hard_ += msg->s2_hard<maxint_t>();

  if (msg->low() > max_low_)
  {
    max_low_ = msg->low();
    segments_ = msg->segments();
    segment_size_ = msg->segment_size();

    Runtime runtime;
    runtime.init = msg->init_seconds();
    runtime.secs = msg->seconds();

    double next = get_next(runtime);
    next = in_between(0.25, next, 2.0);
    next *= segments_;
    next = max(1.0, next);
    segments_ = (int64_t) next;
  }

  auto high = low_ + segments_ * segment_size_;

  // Most hard special leaves are located just past
  // smallest_hard_leaf_. In order to prevent assigning
  // the bulk of work to a single thread we reduce
  // the number of segments to a minimum.
  //
  if (smallest_hard_leaf_ >= low_ &&
      smallest_hard_leaf_ <= high)
  {
    segments_ = 1;
  }

  // udpate msg with new work todo
  msg->update(low_, segments_, segment_size_);

  low_ += segments_ * segment_size_;
  low_ = min(low_, limit_);
}

double MpiLoadBalancer::get_next(Runtime& runtime) const
{
  double min_secs = runtime.init * 10;
  double run_secs = runtime.secs;

  min_secs = max(min_secs, 0.01);
  run_secs = max(run_secs, min_secs / 10);

  double rem = remaining_secs();
  double threshold = rem / 4;
  threshold = max(threshold, min_secs);

  return threshold / run_secs;
}

/// Remaining seconds till finished
double MpiLoadBalancer::remaining_secs() const
{
  double percent = status_.getPercent(low_, z_, s2_hard_, s2_approx_);
  percent = in_between(20, percent, 100);

  double total_secs = get_wtime() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

} // namespace
