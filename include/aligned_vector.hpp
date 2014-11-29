///
/// @file  aligned_vector.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
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
  #define CACHE_LINE_SIZE 128
#endif

namespace primecount {

/// The aligned_vector class aligns each of its elements on a
/// new cache line in order to avoid false sharing (cache trashing)
/// when multiple threads write to adjacent elements.
///
template <typename T, std::size_t ALIGN = CACHE_LINE_SIZE>
class aligned_vector
{
public:
  aligned_vector(std::size_t size)
    : vector_(size) { }
  T& operator[](std::size_t pos) { return vector_[pos].val; }
  std::size_t size() const { return vector_.size(); }
private:
  struct aligned_type
  {
    T val;
    char pad[ALIGN - sizeof(T)];
  };
  std::vector<aligned_type> vector_;
};

} // namespace

#endif
