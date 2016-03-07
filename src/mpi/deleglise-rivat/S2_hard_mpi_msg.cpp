///
/// @file   S2_hard_mpi_msg.cpp
/// @brief  MPI utility functions and classes.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <S2_hard_mpi_msg.hpp>

#include <cstddef>
#include <mpi.h>
#include <int128.hpp>

#ifndef MPI_INT64_T
  #define MPI_INT64_T MPI_LONG_LONG
#endif

namespace primecount {

S2_hard_mpi_msg::S2_hard_mpi_msg()
{
  init_MPI_struct();
  reset();
}

S2_hard_mpi_msg::S2_hard_mpi_msg(int64_t proc_id,
                                 int64_t low,
                                 int64_t high,
                                 int64_t segment_size,
                                 int64_t segments_per_thread)
{
  init_MPI_struct();
  set(proc_id, low, high, segment_size, segments_per_thread);
}

void S2_hard_mpi_msg::set(int64_t proc_id,
                          int64_t low,
                          int64_t high,
                          int64_t segment_size,
                          int64_t segments_per_thread)
{
  reset();
  msgData_.proc_id = proc_id;
  msgData_.low = low;
  msgData_.high = high;
  msgData_.segment_size = segment_size;
  msgData_.segments_per_thread = segments_per_thread;
}

void S2_hard_mpi_msg::reset()
{
  msgData_.proc_id = 0;
  msgData_.low = 0;
  msgData_.high = 0;
  msgData_.segment_size = 0;
  msgData_.segments_per_thread = 0;
  msgData_.s2_hard[0] = 0;
  msgData_.s2_hard[1] = 0;
  msgData_.init_seconds = 0;
  msgData_.seconds = 0;
  msgData_.rsd = 0;
  msgData_.finished = false;
}

S2_hard_mpi_msg::~S2_hard_mpi_msg()
{
  MPI_Type_free(&mpi_type_);
}

void S2_hard_mpi_msg::init_MPI_struct()
{
  int items = 10;
  int block_lengths[10] = { 1, 1, 1, 1, 1, 2, 1, 1, 1, 1 };

  MPI_Datatype types[10] = { MPI_INT,     // proc_id
                            MPI_INT64_T, // low
                            MPI_INT64_T, // high
                            MPI_INT64_T, // segment_size
                            MPI_INT64_T, // segments_per_thread
                            MPI_INT64_T, // s2_hard
                            MPI_DOUBLE,  // init_seconds
                            MPI_DOUBLE,  // seconds
                            MPI_DOUBLE,  // rsd
                            MPI_INT };   // finished
  MPI_Aint offsets[10];

  offsets[0] = offsetof(MsgData, proc_id);
  offsets[1] = offsetof(MsgData, low);
  offsets[2] = offsetof(MsgData, high);
  offsets[3] = offsetof(MsgData, segment_size);
  offsets[4] = offsetof(MsgData, segments_per_thread);
  offsets[5] = offsetof(MsgData, s2_hard);
  offsets[6] = offsetof(MsgData, init_seconds);
  offsets[7] = offsetof(MsgData, seconds);
  offsets[8] = offsetof(MsgData, rsd);
  offsets[9] = offsetof(MsgData, finished);

  MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_type_);
  MPI_Type_commit(&mpi_type_);
}

void S2_hard_mpi_msg::send(int proc_id)
{
  MPI_Send(&msgData_, 1, mpi_type_, proc_id, proc_id, MPI_COMM_WORLD);
}

void S2_hard_mpi_msg::send_finish()
{
  msgData_.finished = true;
  int proc_id = msgData_.proc_id;

  MPI_Send(&msgData_, 1, mpi_type_, proc_id, proc_id, MPI_COMM_WORLD);
}

void S2_hard_mpi_msg::recv(int proc_id)
{
  MPI_Status status;
  MPI_Recv(&msgData_, 1, mpi_type_, mpi_master_proc_id(), proc_id, MPI_COMM_WORLD, &status);
}

void S2_hard_mpi_msg::recv_any()
{
  MPI_Status status;
  MPI_Recv(&msgData_, 1, mpi_type_, MPI_ANY_SOURCE, mpi_master_proc_id(), MPI_COMM_WORLD, &status);
}

int S2_hard_mpi_msg::proc_id() const
{
  return msgData_.proc_id;
}

int64_t S2_hard_mpi_msg::low() const
{
  return msgData_.low;
}

int64_t S2_hard_mpi_msg::high() const
{
  return msgData_.high;
}

int64_t S2_hard_mpi_msg::segment_size() const
{
  return msgData_.segment_size;
}

int64_t S2_hard_mpi_msg::segments_per_thread() const
{
  return msgData_.segments_per_thread;
}

double S2_hard_mpi_msg::init_seconds() const
{
  return msgData_.init_seconds;
}

double S2_hard_mpi_msg::seconds() const
{
  return msgData_.seconds;
}

double S2_hard_mpi_msg::rsd() const
{
  return msgData_.rsd;
}

bool S2_hard_mpi_msg::finished() const
{
  return ((bool) msgData_.finished) == true;
}

} // namespace
