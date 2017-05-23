///
/// @file   MpiMsg.cpp
/// @brief  MpiMsg is used to send and receive
///         messages between the MPI master process and the
///         MPI slave processes during the computation of
///         the hard special leaves.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <MpiMsg.hpp>
#include <int128_t.hpp>

#include <cstddef>
#include <limits>
#include <mpi.h>

#ifndef MPI_INT64_T
  #define MPI_INT64_T MPI_LONG_LONG
#endif

namespace primecount {

MpiMsg::MpiMsg()
{
  init_MPI_struct();

  msgData_.proc_id = 0;
  msgData_.thread_id = 0;
  msgData_.low = 0;
  msgData_.segments = 0;
  msgData_.segment_size = 0;
  msgData_.s2_hard[0] = 0;
  msgData_.s2_hard[1] = 0;
  msgData_.init_seconds = 0;
  msgData_.seconds = 0;
  msgData_.finished = 0;
}

void MpiMsg::set_finished()
{
  msgData_.finished = true;
}

void MpiMsg::update(int64_t low,
                    int64_t segments,
                    int64_t segment_size)
{
  msgData_.low = low;
  msgData_.segments = segments;
  msgData_.segment_size = segment_size;
}

MpiMsg::~MpiMsg()
{
  MPI_Type_free(&mpi_type_);
}

void MpiMsg::init_MPI_struct()
{
  int items = 9;
  int block_lengths[9] = { 1, 1, 1, 1, 1, 2, 1, 1, 1 };

  MPI_Datatype types[9] = { MPI_INT,     // proc_id
                            MPI_INT,     // thread_id
                            MPI_INT64_T, // low
                            MPI_INT64_T, // segments
                            MPI_INT64_T, // segment_size
                            MPI_INT64_T, // s2_hard
                            MPI_DOUBLE,  // init_seconds
                            MPI_DOUBLE,  // seconds
                            MPI_INT };   // finished
  MPI_Aint offsets[9];

  offsets[0] = offsetof(MsgData, proc_id);
  offsets[1] = offsetof(MsgData, thread_id);
  offsets[2] = offsetof(MsgData, low);
  offsets[3] = offsetof(MsgData, segments);
  offsets[4] = offsetof(MsgData, segment_size);
  offsets[5] = offsetof(MsgData, s2_hard);
  offsets[6] = offsetof(MsgData, init_seconds);
  offsets[7] = offsetof(MsgData, seconds);
  offsets[8] = offsetof(MsgData, finished);

  MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_type_);
  MPI_Type_commit(&mpi_type_);
}

void MpiMsg::send(int proc_id)
{
  MPI_Send(&msgData_, 1, mpi_type_, proc_id, proc_id, MPI_COMM_WORLD);
}

void MpiMsg::recv(int proc_id)
{
  MPI_Status status;
  MPI_Recv(&msgData_, 1, mpi_type_, mpi_master_proc_id(), proc_id, MPI_COMM_WORLD, &status);
}

void MpiMsg::recv_any()
{
  MPI_Status status;
  MPI_Recv(&msgData_, 1, mpi_type_, MPI_ANY_SOURCE, mpi_master_proc_id(), MPI_COMM_WORLD, &status);
}

int MpiMsg::proc_id() const
{
  return msgData_.proc_id;
}

int MpiMsg::thread_id() const
{
  return msgData_.thread_id;
}

int MpiMsg::finished() const
{
  return msgData_.finished;
}

int64_t MpiMsg::low() const
{
  return msgData_.low;
}

int64_t MpiMsg::segments() const
{
  return msgData_.segments;
}

int64_t MpiMsg::segment_size() const
{
  return msgData_.segment_size;
}

double MpiMsg::init_seconds() const
{
  return msgData_.init_seconds;
}

double MpiMsg::seconds() const
{
  return msgData_.seconds;
}

} // namespace
