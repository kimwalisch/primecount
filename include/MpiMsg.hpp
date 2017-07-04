///
/// @file  MpiMsg.hpp
/// @brief MpiMsg is used to send and receive messages between
///        the MPI master process and the MPI slave processes
///        during the computation of the hard special leaves.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MPIMSG_HPP
#define MPIMSG_HPP

#include <mpi.h>
#include <int128_t.hpp>
#include <cassert>

namespace primecount {

class MpiMsg
{
public:
  MpiMsg();
  ~MpiMsg();
  void init_MPI_struct();
  void send(int proc_id);
  void recv(int proc_id);
  void recv_any();
  void set_finished();
  int proc_id() const;
  int thread_id() const;
  int finished() const;
  int64_t low() const;
  int64_t segments() const;
  int64_t segment_size() const;
  double init_seconds() const;
  double seconds() const;

  void update(int64_t low,
              int64_t segments,
              int64_t segment_size);

  template <typename T>
  void set(int proc_id,
           int thread_id,
           int64_t low,
           int64_t segments,
           int64_t segment_size,
           T s2_hard,
           double init_seconds,
           double seconds)
  {
    msgData_.proc_id = proc_id;
    msgData_.thread_id = thread_id;
    msgData_.low = low;
    msgData_.segments = segments;
    msgData_.segment_size = segment_size;
    *((T*) msgData_.s2_hard) = s2_hard;
    msgData_.init_seconds = init_seconds;
    msgData_.seconds = seconds;
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
    int thread_id;
    int64_t low;
    int64_t segments;
    int64_t segment_size;
    int64_t s2_hard[2];
    double init_seconds;
    double seconds;
    int finished;
  };

  MsgData msgData_;
  MPI_Datatype mpi_type_;
};

} // namespace

#endif
