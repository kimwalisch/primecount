///
/// @file  PhiCache.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHICACHE_HPP
#define PHICACHE_HPP

#include "PhiTiny.hpp"

#include <stdint.h>
#include <vector>

namespace primecount {

class PhiCache {
public:
  PhiCache(const std::vector<int32_t>& primes, const PhiTiny& phiTiny);
  int64_t phi(int64_t x, int64_t a);
  int64_t phi(int64_t x, int64_t a, int sign);
private:
  /// Cache for phi(x, a) results
  std::vector<std::vector<uint16_t> > cache_;
  const std::vector<int32_t>& primes_;
  const PhiTiny& phiTiny_;
  int64_t bytes_;

  int64_t phi_bsearch(int64_t x, int64_t a) const;
  bool is_phi_bsearch(int64_t x, int64_t a) const;
  bool is_cached(int64_t x, int64_t a) const;
  bool write_to_cache(int64_t x, int64_t a);
};

} // namespace primecount

#endif
