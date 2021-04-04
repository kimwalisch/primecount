///
/// @file  StatusAC.cpp
/// @brief The StatusAC class is used to print the status (in percent)
///        of the A & C formulas in Xavier Gourdon's algorithm.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <StatusAC.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <iostream>

#if defined(_OPENMP)
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace primecount {

bool StatusAC::isPrint(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

void StatusAC::print(int64_t low, int64_t limit, int64_t segment_size)
{
  // check --status option used
  if (!is_print())
    return;

#if defined(_OPENMP)
  if (omp_get_thread_num() != 0)
    return;
#endif

  double time = get_time();

  if (isPrint(time))
  {
    time_ = time;
    int64_t high = low + segment_size;
    cout << "\rSegment: " << high / segment_size << "/" << ceil_div(limit, segment_size) << flush;
  }
}

} // namespace
