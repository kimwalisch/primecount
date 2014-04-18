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

class PhiCache
{
public:
  PhiCache(const std::vector<int32_t>& primes, const PhiTiny& phiTiny);
  int64_t phi(int64_t x, int64_t a);
  int64_t phi(int64_t x, int64_t a, int sign);
private:
  enum {
    /// Cache phi(x, a) results if a <= CACHE_A_LIMIT
    CACHE_A_LIMIT = 500,
    /// Keep the cache size below CACHE_BYTES_LIMIT per thread
    CACHE_BYTES_LIMIT = 16 << 20
  };
  /// Cache of phi(x, a) results
  std::vector<std::vector<uint16_t> > cache_;
  const std::vector<int32_t>& primes_;
  const PhiTiny& phiTiny_;
  int64_t bytes_;

  /// Disable the default assignment operator
  void operator=(const PhiCache&);
  int64_t phi_bsearch(int64_t x, int64_t a) const;
  bool is_phi_bsearch(int64_t x, int64_t a) const;
  bool write_to_cache(int64_t x, int64_t a);

  bool is_cached(int64_t x, int64_t a) const
  {
    return a <= CACHE_A_LIMIT &&
           x < static_cast<int64_t>(cache_[a].size()) &&
           cache_[a][x] != 0;
  }
};

} // namespace primecount

#endif
