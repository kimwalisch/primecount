///
/// @file   S2_hard_mpi_LoadBalancer.cpp
/// @brief  The S2_hard_mpi_LoadBalancer evenly distributes the
///         computation of the hard special leaves onto cluster nodes.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2_hard_mpi_LoadBalancer.hpp>
#include <S2_hard_mpi_msg.hpp>
#include <primecount-internal.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;

namespace primecount {

S2_hard_mpi_LoadBalancer::S2_hard_mpi_LoadBalancer(int64_t high,
                                                   int64_t y,
                                                   int64_t z,
                                                   int slave_procs) :
  low_(0),
  high_(high),
  y_(y),
  z_(z),
  slave_procs_(slave_procs),
  max_finished_(0),
  segment_size_(isqrt(z)),
  segments_per_thread_(1),
  proc_interval_(0),
  rsd_(0),
  start_time_(get_wtime()),
  init_seconds_(0),
  seconds_(0)
{ }

bool S2_hard_mpi_LoadBalancer::finished() const
{
  return low_ > z_;
}

bool S2_hard_mpi_LoadBalancer::is_increase(double percent) const
{
  double min_seconds = max(0.1, init_seconds_ * 10);

  if (seconds_ < min_seconds)
    return true;

  // avoid division by 0
  percent = in_between(1, percent, 100);

  // calculate remaining time till finished
  double elapsed_time = get_wtime() - start_time_;
  double remaining_time = elapsed_time * (100 / percent) - elapsed_time;
  double max_seconds = remaining_time / 4;
  double is_increase = max(min_seconds, max_seconds);

  return seconds_ < is_increase;
}

void S2_hard_mpi_LoadBalancer::update(S2_hard_mpi_msg* msg,
                                      double percent)
{
  if (msg->high() >= max_finished_)
  {
    max_finished_ = msg->high();
    proc_interval_ = msg->high() - msg->low();
    segment_size_ = msg->segment_size();
    segments_per_thread_ = msg->segments_per_thread();
    rsd_ = msg->rsd();
    init_seconds_ = msg->init_seconds();
    seconds_ = msg->seconds();
  }

  int64_t next_interval = proc_interval_;

  // balance load by increasing or decreasing the next
  // interval based on previous run-time
  if (is_increase(percent))
    next_interval *= 2;
  else
    next_interval = max(next_interval / 2, isqrt(z_)); 

  low_ = high_ + 1;
  high_ = min(low_ + next_interval, z_);

  // udpate existing message with new work todo
  msg->set(msg->proc_id(), low_, high_, segment_size_, segments_per_thread_, rsd_);
}

} // namespace
