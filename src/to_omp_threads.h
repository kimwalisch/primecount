///
/// @file  to_omp_threads.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef TO_OMP_THREADS_H
#define TO_OMP_THREADS_H

#ifdef _OPENMP

#include <primecount.h>
#include <omp.h>

namespace primecount {

inline int to_omp_threads(int threads)
{
  return (threads != MAX_THREADS) ? threads : omp_get_max_threads();
}

} // namespace primecount

#endif /* _OPENMP */

#endif
