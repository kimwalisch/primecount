///
/// @file  aligned_vector.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ALIGNED_VECTOR_HPP
#define ALIGNED_VECTOR_HPP

#include <pod_vector.hpp>
#include <cstddef>

/// Maximum cache line size of current CPUs.
/// Note that in 2019 all x86 CPU have a cache line size of 64 bytes.
/// However there are CPUs out there that have much larger cache line
/// sizes e.g. IBM z13 CPUs from 2015 have a cache line size of 256
/// bytes. Hence in order to be future-proof we set the maximum cache
/// line size to 512 bytes.
///
#ifndef CACHE_LINE_SIZE
  #define CACHE_LINE_SIZE (1 << 9)
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
  static_assert(sizeof(T) < CACHE_LINE_SIZE,
                "sizeof(T) must be < CACHE_LINE_SIZE");

public:
  void resize(std::size_t size)
  {
    vect_.resize(size);

    // Default initialize values.
    // Do not initialize the padding memory.
    for (std::size_t i = 0; i < size; i++)
      vect_[i].val = {};
  }
  char unused()
  {
    vect_[0].pad[0] = 123;
    return vect_[0].pad[0];
  }

  aligned_vector() = default;
  aligned_vector(std::size_t size) { resize(size); }
  std::size_t size() const { return vect_.size(); }
  T& operator[](std::size_t pos) { return vect_[pos].val; }

private:
  struct CacheLine {
    T val;
    // We cannot use alignas(CACHE_LINE_SIZE) for the
    // CacheLine struct as GCC does not yet support
    // alignas(n) with n > 128. Also alignas(n) for
    // over-aligned data and dynamic memory allocation
    // is only supported since C++17.
    char pad[CACHE_LINE_SIZE - sizeof(T)];
  };

  pod_vector<CacheLine> vect_;
};

} // namespace

#endif
