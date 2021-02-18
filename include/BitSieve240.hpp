///
/// @file  BitSieve240.hpp
/// @brief The BitSieve240 base class contains lookup tables that are
///        needed to implement a prime sieving algorithm where each
///        bit corresponds to an integer that is not divisible by 2,
///        3 and 5. The 8 bits of each byte represent an interval of
///        size 30. Hence one uint64_t sieve array element (8 bytes)
///        corresponds to an interval of 30 * 8 = 240.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE240_HPP
#define BITSIEVE240_HPP

#include <stdint.h>
#include <array>

namespace primecount {

class BitSieve240
{
protected:
  static const std::array<uint64_t, 240> unset_bit_;
  static const std::array<uint64_t, 240> unset_larger_;
};

} // namespace

#endif
