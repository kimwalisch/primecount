///
/// @file  S2.hpp.
/// @brief S2 function declarations.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HPP
#define S2_HPP

#include <int128.hpp>
#include <stdint.h>

namespace primecount {

/// ------------------------ S2_trivial() ----------------------------

int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   int threads);

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads);

#endif

/// ------------------------ S2_easy() -------------------------------

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int threads);

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads);

#endif

#ifdef HAVE_MPI

int64_t S2_easy_mpi(int64_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads);

#ifdef HAVE_INT128_T

int128_t S2_easy_mpi(int128_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     int threads);

#endif

#endif

/// ------------------------ S2_hard() -------------------------------

int64_t S2_hard(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int64_t s2_hard_approx,
                int threads);

#ifdef HAVE_INT128_T

int128_t S2_hard(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int128_t s2_hard_approx,
                 int threads);

#endif

#ifdef HAVE_MPI

int64_t S2_hard_mpi(int64_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int64_t s2_hard_approx,
                    int threads);

#ifdef HAVE_INT128_T

int128_t S2_hard_mpi(int128_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     int128_t s2_hard_approx,
                     int threads);

#endif

#endif

} // namespace

#endif
