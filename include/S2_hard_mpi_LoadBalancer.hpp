///
/// @file   S2_hard_mpi_LoadBalancer.hpp
/// @brief  The S2_hard_mpi_LoadBalancer evenly distributes
///         the computation of the hard special leaves onto
///         cluster nodes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HARD_MPI_LOADBALANCER_HPP
#define S2_HARD_MPI_LOADBALANCER_HPP

#include <mpi.h>
#include <S2_hard_mpi_msg.hpp>
#include <S2Status.hpp>
#include <int128_t.hpp>
#include <stdint.h>

namespace primecount {

class S2_hard_mpi_LoadBalancer
{
public:
  S2_hard_mpi_LoadBalancer(maxint_t x, int64_t z, int64_t high, maxint_t s2_approx);
  void update(S2_hard_mpi_msg* msg, maxint_t s2_hard);
  bool finished() const;

private:
  bool is_increase(maxint_t s2_hard) const;
  double remaining_secs(maxint_t s2_hard) const;

  int64_t z_;
  int64_t low_;
  int64_t high_;
  int64_t max_finished_;
  int64_t segment_size_;
  int64_t segments_per_thread_;
  int64_t proc_distance_;
  maxint_t s2_approx_;
  double rsd_;
  double start_time_;
  double init_seconds_;
  double seconds_;
  S2Status status_;
};

} // namespace

#endif
