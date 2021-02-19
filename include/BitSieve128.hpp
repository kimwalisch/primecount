///
/// @file  BitSieve128.hpp
/// @brief The BitSieve128 base class contains lookup tables that are
///        needed to implement a prime sieving algorithm where each bit
///        corresponds to an odd integer. The BitSieve128 class uses
///        the uint64_t data type for its sieve array, hence each sieve
///        array element corresponds to an interval of 64 * 2 = 128.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE128_HPP
#define BITSIEVE128_HPP

#include <stdint.h>
#include <array>

namespace primecount {

class BitSieve128
{
protected:
  static const std::array<uint64_t, 128> unset_bit_;
  static const std::array<uint64_t, 128> unset_larger_;
};

} // namespace

#endif
