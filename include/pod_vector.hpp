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

#include <vector>

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
  pod_vector() = default;
  pod_vector(std::size_t size) : vect_(size) { }
  void resize(std::size_t size) { vect_.resize(size); }
  std::size_t size() const { return vect_.size(); }
  std::size_t capacity() const { return vect_.capacity(); }
  T* data() { return (T*) vect_.data(); }
  const T* data() const { return (T*) vect_.data(); }

  // Some Linux distributions (e.g. Fedora) compile primecount
  // with -D_GLIBCXX_ASSERTIONS which enables std::vector
  // bounds checks when using the operator[]. This deteriorates
  // primecount's performance by up to 10%. Therefore, we avoid
  // using std::vector's operator[].
  T& operator[](std::size_t pos) { return (vect_.data() + pos)->val; }
  const T& operator[](std::size_t pos) const { return (vect_.data() + pos)->val; }

private:
  struct NoInitType
  {
    NoInitType() { };
    T val;
  };
  std::vector<NoInitType> vect_;
};

} // namespace

#endif
