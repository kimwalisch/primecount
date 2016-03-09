///
/// @file   S2_hard_mpi_msg.hpp
/// @brief  The S2_hard MPI master process uses the S2_hard_mpi_msg
///         class to send messages to the S2_hard MPI slave processes
///         and assign them work to do. The S2_hard MPI slave
///         processes use the S2_hard_mpi_msg class to send back their
///         partial results to the S2_hard MPI master process.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HARD_MPI_MSG_HPP
#define S2_HARD_MPI_MSG_HPP

#include <mpi.h>
#include <int128.hpp>
#include <cassert>

namespace primecount {

class S2_hard_mpi_msg
{
public:
  S2_hard_mpi_msg();
  S2_hard_mpi_msg(int64_t proc_id,
                  int64_t low,
                  int64_t high,
                  int64_t segment_size,
                  int64_t segments_per_thread);

  ~S2_hard_mpi_msg();
  void reset();
  void init_MPI_struct();
  void send(int proc_id);
  void send_finish();
  void recv(int proc_id);
  void recv_any();
  int proc_id() const;
  int64_t low() const;
  int64_t high() const;
  int64_t segment_size() const;
  int64_t segments_per_thread() const;
  double init_seconds() const;
  double seconds() const;
  double rsd() const;
  bool finished() const;

  void set(int64_t proc_id,
           int64_t low,
           int64_t high,
           int64_t segment_size,
           int64_t segments_per_thread,
           double rsd);

  template <typename T>
  S2_hard_mpi_msg(int proc_id,
                  int64_t low,
                  int64_t high,
                  int64_t segment_size,
                  int64_t segments_per_thread,
                  T s2_hard,
                  double init_seconds,
                  double seconds,
                  double rsd)
  {
    init_MPI_struct();

    msgData_.proc_id = proc_id;
    msgData_.low = low;
    msgData_.high = high;
    msgData_.segment_size = segment_size;
    msgData_.segments_per_thread = segments_per_thread;
    *((T*) msgData_.s2_hard) = s2_hard;
    msgData_.init_seconds = init_seconds;
    msgData_.seconds = seconds;
    msgData_.rsd = rsd;
    msgData_.finished = false;
  }

  template <typename T>
  T s2_hard() const
  {
    return *((T*) msgData_.s2_hard);
  }

private:
  struct MsgData
  {
    int proc_id;
    int64_t low;
    int64_t high;
    int64_t segment_size;
    int64_t segments_per_thread;
    int64_t s2_hard[2];
    double init_seconds;
    double seconds;
    double rsd;
    int finished;
  };

  MsgData msgData_;
  MPI_Datatype mpi_type_;
};

} // namespace

#endif /* S2_HARD_MPI_MSG_HPP */
