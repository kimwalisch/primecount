///
/// @file  validate_threads.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef VALIDATE_THREADS_HPP
#define VALIDATE_THREADS_HPP

#include <primecount.hpp>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

inline int validate_threads(int threads)
{
#ifdef _OPENMP
  if (threads == primecount::MAX_THREADS) 
    threads = omp_get_max_threads();
#endif
  return std::max(1, threads);
}

#endif
