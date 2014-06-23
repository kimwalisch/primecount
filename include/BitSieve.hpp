///
/// @file  BitSieve.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE_HPP
#define BITSIEVE_HPP

#include <primecount-internal.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>
#include <stdint.h>

const unsigned unset_bit[32] =
{
  ~(1u <<  0), ~(1u <<  1), ~(1u <<  2),
  ~(1u <<  3), ~(1u <<  4), ~(1u <<  5),
  ~(1u <<  6), ~(1u <<  7), ~(1u <<  8),
  ~(1u <<  9), ~(1u << 10), ~(1u << 11),
  ~(1u << 12), ~(1u << 13), ~(1u << 14),
  ~(1u << 15), ~(1u << 16), ~(1u << 17),
  ~(1u << 18), ~(1u << 19), ~(1u << 20),
  ~(1u << 21), ~(1u << 22), ~(1u << 23),
  ~(1u << 24), ~(1u << 25), ~(1u << 26),
  ~(1u << 27), ~(1u << 28), ~(1u << 29),
  ~(1u << 30), ~(1u << 31)
};

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

  /// Count the number of 1 bits inside [start, stop]
  int64_t count(int64_t start, int64_t stop) const
  {
    assert(stop < (int64_t) size_);
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(bits_.data());
    return popcount(bits, start, stop);
  }

  /// Set all bits to 1, except bits corresponding
  /// to 0, 1 and even numbers > 2.
  /// 
  void memset(uint64_t low)
  {
    unsigned mask = (low & 1) ? 0x55555555u : 0xAAAAAAAAu;
    std::fill(bits_.begin(), bits_.end(), mask);

    // correct 0, 1 and 2
    if (low <= 2)
    {
      uint32_t bitmask = ~((1u << (2 - low)) - 1);
      bits_[0] &= bitmask;
      bits_[0] |= 1 << (2 - low);
    }
  }

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
    bits_[pos >> 5] &= unset_bit[pos & 31];
  }
private:
  std::vector<uint32_t> bits_;
  std::size_t size_;
};

} // namespace

#endif
