///
/// @file  PhiTiny.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <stdint.h>
#include <vector>

namespace primecount {

class PhiTiny {
public:
  PhiTiny();
  int64_t phi(int64_t x, int64_t a) const;
  static int64_t max_a() { return 6; }
private:
  std::vector<int16_t> phi_cache_[7];
};

inline bool is_phi_tiny(int64_t a)
{
  return a <= PhiTiny::max_a();
}

int64_t phi_tiny(int64_t x, int64_t a);

} // namespace primecount

#endif
