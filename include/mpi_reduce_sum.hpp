///
/// @file  mpi_reduce_sum.hpp
/// @brief MPI utility functions and classes
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MPI_REDUCE_SUM_HPP
#define MPI_REDUCE_SUM_HPP

#include <mpi.h>
#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <cassert>

#ifndef MPI_INT64_T
  #define MPI_INT64_T MPI_LONG_LONG
#endif

namespace primecount {

template <typename T>
inline void mpi_sum_helper(int* x, int* sum)
{
  T& x_ = *((T*) x);
  T& sum_ = *((T*) sum);
  sum_ += x_;
}

inline void mpi_sum(int* x, int* sum, int* len, MPI_Datatype *dtype)
{
  unused_param(dtype);
  if (*len * sizeof(int64_t) == sizeof(int64_t))
    mpi_sum_helper<int64_t>(x, sum);
#if defined(HAVE_INT128_T)
  if (*len * sizeof(int64_t) == sizeof(int128_t))
    mpi_sum_helper<int128_t>(x, sum);
#endif
}

/// Works with both int64_t and int128_t
template <typename T>
inline T mpi_reduce_sum(T x)
{
  assert(sizeof(T) >= sizeof(int64_t));
  T sum = 0;

  MPI_Op op;
  MPI_Op_create((MPI_User_function*) mpi_sum, 1, &op);
  MPI_Reduce((int*) &x, (int*) &sum, sizeof(T) / sizeof(int64_t), MPI_INT64_T, op, 0, MPI_COMM_WORLD);
  MPI_Op_free(&op);

  return sum;
}

} // namespace

#endif
