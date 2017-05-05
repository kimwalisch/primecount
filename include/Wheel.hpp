///
/// @file   Wheel.hpp
/// @brief  Data structures related to wheel factorization.
///         Wheel factorization is used to skip multiples of
///         small primes in the sieve of Eratosthenes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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

/// The InitWheel data structure is used to calculate the
/// first multiple >= start of a sieving prime
///
struct InitWheel
{
  int8_t next_multiple_factor;
  int8_t wheel_index;
};

/// The NextWheel data structure is used to calculate the
/// next multiple of a sieving prime
///
struct NextWheel
{
  int8_t next_multiple_factor;
  int8_t next_wheel_index;
};

/// For each sieving prime we create a WheelItem which
/// contains the sieving prime's next multiple and the
/// wheel index of that multiple
///
struct WheelItem
{
  WheelItem(int64_t multiple, int64_t index) :
    next_multiple(multiple),
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

/// 4th wheel, skips multiples of 2, 3, 5 and 7
class Wheel
{
public:
  /// Calculate the first multiple >= low of each prime.
  /// When sieving special leaves both multiples
  /// and primes are crossed-off.
  ///
  template <typename Primes>
  Wheel(Primes& primes,
        int64_t size,
        int64_t low)
  {
    wheel_.reserve(size);
    wheel_.emplace_back(0, 0);

    for (int64_t b = 1; b < size; b++)
    {
      int64_t prime = primes[b];
      int64_t quotient = ceil_div(low, prime);

      // first multiple >= low
      int64_t multiple = prime * quotient;

      // calculate the next multiple of prime that
      // is not divisible by any of the wheel's
      // prime factors (2, 3, 5, 7)
      int64_t next_multiple_factor = init[quotient % 210].next_multiple_factor;
      int64_t wheel_index = init[quotient % 210].wheel_index;
      multiple += prime * next_multiple_factor;

      wheel_.emplace_back(multiple, wheel_index);
    }
  }

  static int64_t next_multiple_factor(int64_t* wheel_index)
  {
    int64_t next_multiple_factor = next[*wheel_index].next_multiple_factor;
    *wheel_index = next[*wheel_index].next_wheel_index;
    return next_multiple_factor;
  }

  WheelItem& operator[](int64_t i)
  {
    return wheel_[i];
  }
private:

  static const InitWheel init[210];
  static const NextWheel next[48];
  std::vector<WheelItem> wheel_;
};

} // namespace

#endif
