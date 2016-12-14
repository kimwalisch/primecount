///
/// @file   Wheel.hpp
/// @brief  Data structures related to wheel factorization.
///         Wheel factorization is used to skip multiples of small
///         primes in the sieve of Eratosthenes.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef WHEEL_HPP
#define WHEEL_HPP

#include <imath.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

/// The InitWheel data structure is used to calculate the first
/// multiple >= start of each sieving prime.
///
struct InitWheel
{
  int8_t next_multiple_factor;
  int8_t wheel_index;
};

struct NextWheel
{
  int8_t next_multiple_factor;
  int8_t next_wheel_index;
};

struct WheelItem
{
  WheelItem(int64_t multiple, int64_t index)
    : next_multiple(multiple), 
      wheel_index((int8_t) index)
  { }

  void set(int64_t multiple, 
           int64_t next_wheel_index)
  {
    next_multiple = multiple;
    wheel_index = (int8_t) next_wheel_index;
  }

  int64_t next_multiple;
  int8_t wheel_index;
};

/// 4th wheel, skips multiples of 2, 3, 5 and 7.
class Wheel
{
public:
  /// Calculate the first multiple >= low of each prime.
  /// When sieving special leaves both multiples and
  /// primes are crossed-off.
  ///
  template <typename Primes>
  Wheel(Primes& primes,
        int64_t size,
        int64_t low)
  {
    wheelItems_.reserve(size);
    push_back(0, 0);

    for (int64_t b = 1; b < size; b++)
    {
      int64_t prime = primes[b];
      int64_t quotient = ceil_div(low, prime);

      // calculate the first multiple of prime >= low
      int64_t multiple = prime * quotient;

      // calculate the next multiple of prime that is not
      // divisible by any of the wheel's factors (2, 3, 5, 7)
      int64_t next_multiple_factor = initWheel210[quotient % 210].next_multiple_factor;
      int64_t wheel_index = initWheel210[quotient % 210].wheel_index;
      multiple += prime * next_multiple_factor;

      push_back(multiple, wheel_index);
    }
  }

  /// Calculate the next multiple of prime using:
  /// next_multiple = multiple + prime * next_multiple_factor(&wheel_index)
  ///
  static int64_t next_multiple_factor(int64_t* wheel_index)
  {
    int64_t next_multiple_factor = nextWheel210[*wheel_index].next_multiple_factor;
    *wheel_index = nextWheel210[*wheel_index].next_wheel_index;
    return next_multiple_factor;
  }

  WheelItem& operator[](int64_t i)
  {
    return wheelItems_[i];
  }
private:
  void push_back(int64_t multiple,
                 int64_t wheel_index)
  {
    #if __cplusplus >= 201103L
      wheelItems_.emplace_back(multiple, wheel_index);
    #else
      wheelItems_.push_back(WheelItem(multiple, wheel_index));
    #endif
  }

  static const InitWheel initWheel210[210];
  static const NextWheel nextWheel210[48];
  std::vector<WheelItem> wheelItems_;
};

} // namespace

#endif
