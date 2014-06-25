///
/// @file  BitSieve.hpp
/// @brief Bit array for prime sieving.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE_HPP
#define BITSIEVE_HPP

#include <cassert>
#include <cstddef>
#include <vector>
#include <stdint.h>

namespace primecount {

/// The BitSieve data structure uses bit packing to save memory.
/// BitSieve uses 1 byte of memory for 8 numbers, each bit
/// corresponds to one integer.
///
class BitSieve
{
public:
  BitSieve(std::size_t size)
    : bits_((size + 31) / 32 + sizeof(uint64_t)),
      size_(size)
  { }

  /// Count the number of 1 bits inside the interval [start, stop]
  uint64_t count(uint64_t start, uint64_t stop) const;

  /// Set all bits to 1, except bits corresponding
  /// to 0, 1 and even numbers > 2.
  void memset(uint64_t low);

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    unsigned mask = 1u << (pos & 31);
    return (bits_[pos >> 5] & mask) != 0;
  }

  std::size_t size() const
  {
    return size_;
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    bits_[pos >> 5] &= unset_bit_[pos & 31];
  }
private:
  static const unsigned int unset_bit_[32];
  std::vector<uint32_t> bits_;
  std::size_t size_;
  static uint64_t popcount(const uint64_t*, uint64_t, uint64_t);
  static uint64_t popcount_edges(const uint64_t*, uint64_t, uint64_t);
};

} // namespace

#endif
