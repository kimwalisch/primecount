///
/// @file   mpi_reduce_sum.cpp
/// @brief  Add 128-bit MPI sum reduction support.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <mpi_reduce_sum.hpp>

#include <cstddef>
#include <mpi.h>
#include <int128.hpp>

namespace {

template <typename T>
inline void mpi_sum_t(int* x, int* sum)
{
  T& x_ = *((T*) x);
  T& sum_ = *((T*) sum);
  sum_ += x_;
}

} // namespace

namespace primecount {

void mpi_sum(int* x, int* sum, int* len, MPI_Datatype *dtype)
{
  if (*len * sizeof(int64_t) == sizeof(int64_t))
    mpi_sum_t<int64_t>(x, sum);
#if defined(HAVE_INT128_T)
  if (*len * sizeof(int64_t) == sizeof(int128_t))
    mpi_sum_t<int128_t>(x, sum);
#endif
}

} // namespace
