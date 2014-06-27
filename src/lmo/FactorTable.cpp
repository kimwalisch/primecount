///
/// @file  FactorTable.cpp
/// @see   FactorTable.hpp for documentation.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <FactorTable.hpp>
#include <pmath.hpp>

#include <limits>
#include <stdexcept>
#include <stdint.h>
#include <vector>

using namespace std;

namespace primecount {

const uint8_t FactorTable::numbers_[48] =
{
    1,  11,  13,  17,  19,  23,
   29,  31,  37,  41,  43,  47,
   53,  59,  61,  67,  71,  73,
   79,  83,  89,  97, 101, 103,
  107, 109, 113, 121, 127, 131,
  137, 139, 143, 149, 151, 157,
  163, 167, 169, 173, 179, 181,
  187, 191, 193, 197, 199, 209
};

const int8_t FactorTable::indexes_[210] =
{
  -1,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  1,  1,  2,  2,  2,
   2,  3,  3,  4,  4,  4,  4,  5,
   5,  5,  5,  5,  5,  6,  6,  7,
   7,  7,  7,  7,  7,  8,  8,  8,
   8,  9,  9, 10, 10, 10, 10, 11,
  11, 11, 11, 11, 11, 12, 12, 12,
  12, 12, 12, 13, 13, 14, 14, 14,
  14, 14, 14, 15, 15, 15, 15, 16,
  16, 17, 17, 17, 17, 17, 17, 18,
  18, 18, 18, 19, 19, 19, 19, 19,
  19, 20, 20, 20, 20, 20, 20, 20,
  20, 21, 21, 21, 21, 22, 22, 23,
  23, 23, 23, 24, 24, 25, 25, 25,
  25, 26, 26, 26, 26, 26, 26, 26,
  26, 27, 27, 27, 27, 27, 27, 28,
  28, 28, 28, 29, 29, 29, 29, 29,
  29, 30, 30, 31, 31, 31, 31, 32,
  32, 32, 32, 32, 32, 33, 33, 34,
  34, 34, 34, 34, 34, 35, 35, 35,
  35, 35, 35, 36, 36, 36, 36, 37,
  37, 38, 38, 38, 38, 39, 39, 39,
  39, 39, 39, 40, 40, 41, 41, 41,
  41, 41, 41, 42, 42, 42, 42, 43,
  43, 44, 44, 44, 44, 45, 45, 46,
  46, 46, 46, 46, 46, 46, 46, 46,
  46, 47
};

FactorTable::FactorTable(int64_t max) :
  max_((max < 8) ? 8 : max)
{
  if (isqrt(max_) >= numeric_limits<uint16_t>::max())
    throw runtime_error("FactorTable: sqrt(max) must be < max(uint16_t).");
  init();
}

void FactorTable::init()
{
  factors_.resize(get_index(max_) + 1, 0xffff);
  // mu(1) = 1 -> factors_[0] = lpf - 1
  factors_[0] = 0xffff - 1;

  for (size_t i = 1; i < factors_.size(); i++)
  {
    if (factors_[i] == 0xffff)
    {
      int64_t prime = get_number(i);
      int64_t multiple = prime * get_number(1), j = 2;

      if (prime < 0xffff)
        factors_[i] = (uint16_t) prime;

      for (; multiple <= max_; multiple = prime * get_number(j++))
      {
        int64_t index = get_index(multiple);

        // prime is the smallest factor
        if (factors_[index] == 0xffff)
          factors_[index] = (uint16_t) prime;
        // the least significant bit indicates whether multiple has
        // an even (0) or odd (1) number of prime factors
        else if (factors_[index] > 0)
          factors_[index] ^= 1;
      }

      // Moebius function is 0 if n has a squared prime factor
      multiple = prime * prime * get_number(0), j = 1;
      for (; multiple <= max_; multiple = prime * prime * get_number(j++))
        factors_[get_index(multiple)] = 0;
    }
  }
}

} // namespace
