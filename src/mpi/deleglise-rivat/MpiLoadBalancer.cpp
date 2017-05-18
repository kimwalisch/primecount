///
/// @file   MpiLoadBalancer.cpp
/// @brief  The MpiLoadBalancer evenly distributes the
///         computation of the hard special leaves onto
///         cluster nodes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <MpiLoadBalancer.hpp>
#include <S2_hard_mpi_msg.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;

namespace primecount {

MpiLoadBalancer::MpiLoadBalancer(maxint_t x,
                                 int64_t z,
                                 int64_t high,
                                 maxint_t s2_approx) :
  z_(z),
  low_(0),
  high_(high),
  max_finished_(0),
  segment_size_(isqrt(z)),
  segments_per_thread_(1),
  proc_distance_(0),
  s2_approx_(s2_approx),
  rsd_(0),
  start_time_(get_wtime()),
  init_seconds_(0),
  seconds_(0),
  status_(x)
{ }

bool MpiLoadBalancer::finished() const
{
  return low_ > z_;
}

void MpiLoadBalancer::update(S2_hard_mpi_msg* msg, maxint_t s2_hard)
{
  if (msg->high() >= max_finished_)
  {
    max_finished_ = msg->high();
    proc_distance_ = msg->high() - msg->low();
    segment_size_ = msg->segment_size();
    segments_per_thread_ = msg->segments_per_thread();
    rsd_ = msg->rsd();
    init_seconds_ = msg->init_seconds();
    seconds_ = msg->seconds();
  }

  int64_t distance = proc_distance_;

  // balance load by increasing/decreasing the next
  // distance based on previous run-time
  if (is_increase(s2_hard))
    distance *= 2;
  else
    distance = max(distance / 2, isqrt(z_)); 

  low_ = high_ + 1;
  high_ = min(low_ + distance, z_);

  // udpate existing message with new work todo
  msg->set(msg->proc_id(), low_, high_, segment_size_, segments_per_thread_, rsd_);
}

bool MpiLoadBalancer::is_increase(maxint_t s2_hard) const
{
  double min_secs = max(0.1, init_seconds_ * 10);

  if (seconds_ < min_secs)
    return true;

  double secs = remaining_secs(s2_hard);
  double threshold = secs / 4;
  threshold = max(min_secs, threshold);

  return seconds_ < threshold;
}

/// Remaining seconds till finished
double MpiLoadBalancer::remaining_secs(maxint_t s2_hard) const
{
  double percent = status_.getPercent(low_, z_, s2_hard, s2_approx_);
  percent = in_between(20, percent, 99.9);

  double total_secs = get_wtime() - start_time_;
  double secs = total_secs * (100 / percent) - total_secs;

  return secs;
}

} // namespace
