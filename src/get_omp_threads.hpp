///
/// @file  get_omp_threads.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GET_OMP_THREADS_HPP
#define GET_OMP_THREADS_HPP

#ifdef _OPENMP

#include <primecount.hpp>
#include <algorithm>
#include <omp.h>

namespace primecount {

inline int get_omp_threads(int threads)
{
  return std::max(1, (threads != MAX_THREADS) ? threads : omp_get_max_threads());
}

} // namespace primecount

#endif /* _OPENMP */

#endif
