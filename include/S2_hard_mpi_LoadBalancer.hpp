///
/// @file   S2_hard_mpi_LoadBalancer.hpp
/// @brief  The S2_hard_mpi_LoadBalancer evenly distributes the
///         computation of the hard special leaves onto cluster nodes.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HARD_MPI_LOADBALANCER_HPP
#define S2_HARD_MPI_LOADBALANCER_HPP

#include <mpi.h>
#include <S2_hard_mpi_msg.hpp>
#include <int128.hpp>
#include <cassert>

namespace primecount {

class S2_hard_mpi_LoadBalancer
{
public:
  S2_hard_mpi_LoadBalancer(int64_t low,
                           int64_t high,
                           int64_t y,
                           int64_t z,
                           int slave_procs);

  void update(S2_hard_mpi_msg* msg, double percent);
  bool finished() const;

private:
  bool is_increase(double percent) const;

  int64_t low_;
  int64_t high_;
  int64_t y_;
  int64_t z_;
  int slave_procs_;
  int64_t max_finished_;
  int64_t segment_size_;
  int64_t segments_per_thread_;
  int64_t proc_interval_;
  double start_time_;
  double init_seconds_;
  double seconds_;
};

} // namespace

#endif
