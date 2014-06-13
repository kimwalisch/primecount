///
/// @file   balance_S2_load.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BALANCE_S2_LOAD_HPP
#define BALANCE_S2_LOAD_HPP

#include <aligned_vector.hpp>
#include <stdint.h>

namespace primecount {

void balance_S2_load(int64_t* segment_size,
                     int64_t* segments_per_thread,
                     int64_t min_segment_size,
                     int64_t max_segment_size,
                     double* old_rsd,
                     aligned_vector<double>& timings);

} // namespace primecount

#endif
