///
/// @file  pod_vector.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POD_VECTOR_HPP
#define POD_VECTOR_HPP

#include <primecount.hpp>
#include <cstddef>

namespace primecount {

/// Plain old data vector, does not default initialize memory.
/// Since primecount may allocate gigabytes of memory and
/// afterwards initalize that memory using multiple threads,
/// we don't want our vector to default initialize our memory
/// otherwise we would initialize the same memory twice.
///
template <typename T>
class pod_vector
{
public:
  pod_vector(std::size_t size)
  {
    array_ = new T[size];
    size_ = size;
  }
  /// For performance reasons we don't shrink the memory
  /// when reducing the size of our pod_vector.
  void resize(std::size_t size)
  {
    if (size > size_)
    {
      if (size_ == 0)
        array_ = new T[size];
      else
        throw primecount_error("Increasing the size of non empty pod_vector not supported!");
    }

    size_ = size;
  }
  pod_vector() = default;
  ~pod_vector() { delete[] array_; }
  std::size_t size() const { return size_; }
  T* data() { return array_; }
  const T* data() const { return array_; }
  T& operator[](std::size_t pos) { return array_[pos]; }
  const T& operator[](std::size_t pos) const { return array_[pos]; }

private:
  T* array_ = nullptr;
  size_t size_ = 0;
};

} // namespace

#endif
