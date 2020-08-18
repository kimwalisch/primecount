///
/// @file  MpiMsg.hpp
/// @brief MpiMsg is used to send and receive messages between
///        the MPI main process and the MPI worker processes
///        during the computation of the hard special leaves.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MPIMSG_HPP
#define MPIMSG_HPP

#include <mpi.h>
#include <stdint.h>

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
  maxint_t sum() const;
  double init_seconds() const;
  double seconds() const;

  void update(int64_t low,
              int64_t segments,
              int64_t segment_size);

  void set(int proc_id,
           int thread_id,
           int64_t low,
           int64_t segments,
           int64_t segment_size,
           maxint_t sum,
           double init_seconds,
           double seconds);

private:
  struct MsgData
  {
    int proc_id;
    int thread_id;
    int64_t low;
    int64_t segments;
    int64_t segment_size;
    int64_t sum[2];
    double init_seconds;
    double seconds;
    int finished;
  };

  MsgData msgData_;
  MPI_Datatype mpi_type_;
};

} // namespace

#endif
