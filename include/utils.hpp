///
/// @file  utils.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef UTILS_HPP
#define UTILS_HPP

#include <primecount.hpp>
#include <algorithm>
#include <chrono>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// Get the wall time in seconds.
inline double get_wtime()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  auto now = std::chrono::system_clock::now();
  auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
  return microseconds / 1e6;
#endif
}

inline int validate_threads(int threads)
{
#ifdef _OPENMP
  if (threads == primecount::MAX_THREADS) 
    threads = omp_get_max_threads();
  return std::max(1, threads);
#else
  threads = 1;
  return threads; 
#endif
}

#endif
