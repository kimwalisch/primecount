///
/// @file  BitSieve.hpp
/// @brief The BitSieve class is a bit array for prime sieving.
///        BitSieve packs 8 numbers into one byte i.e. each bit
///        corresponds to one integer.
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

class BitSieve
{
public:
  BitSieve(std::size_t size);

  /// Set all bits to 1, except bits corresponding
  /// to 0, 1 and even numbers > 2.
  void memset(uint64_t low);

  /// Count the number of 1 bits inside the interval [start, stop]
  uint64_t count(uint64_t start, uint64_t stop) const;

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    unsigned bit = 1u << (pos & 31);
    return (bits_[pos >> 5] & bit) != 0;
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    bits_[pos >> 5] &= unset_bit_[pos & 31];
  }

  std::size_t size() const
  {
    return size_;
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
