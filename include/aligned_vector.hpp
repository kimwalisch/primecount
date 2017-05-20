///
/// @file  aligned_vector.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ALIGNED_VECTOR_HPP
#define ALIGNED_VECTOR_HPP

#include <cstddef>
#include <vector>

// Maximum cache line size of current CPUs
#ifndef CACHE_LINE_SIZE
  #define CACHE_LINE_SIZE 512
#endif

namespace primecount {

/// The aligned_vector class aligns each of its
/// elements on a new cache line in order to avoid
/// false sharing (cache trashing) when multiple
/// threads write to adjacent elements
///
template <typename T>
class aligned_vector
{
public:
  aligned_vector(std::size_t size)
    : vect_(size) { }
  std::size_t size() const { return vect_.size(); }
  T& operator[](std::size_t pos) { return vect_[pos].val[0]; }
private:
  struct align_t
  {
    T val[CACHE_LINE_SIZE / sizeof(T)];
  };
  std::vector<align_t> vect_;
};

} // namespace

#endif
