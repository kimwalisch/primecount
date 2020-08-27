///
/// @file  FactorTable.hpp
/// @brief The FactorTable class combines the lpf[n] (least prime
///        factor) and mu[n] (Möbius function) lookup tables into a
///        single factor[n] table which furthermore only contains
///        entries for numbers which are not divisible by 2, 3, 5, 7
///        and 11. The factor[n] lookup table uses up to 19.25
///        times less memory than the lpf[n] & mu[n] lookup tables!
///        factor[n] uses only 2 bytes per entry for 32-bit numbers
///        and 4 bytes per entry for 64-bit numbers.
///
///        The factor table concept was devised and implemented by
///        Christian Bau in 2003. Note that Tomás Oliveira e Silva
///        also suggests combining the mu[n] and lpf[n] lookup tables
///        in his paper. However Christian Bau's FactorTable data
///        structure uses only half as much memory and is also
///        slightly more efficient (uses fewer instructions) than the
///        data structure proposed by Tomás Oliveira e Silva.
///
///        What we store in the factor[n] lookup table:
///
///        1) INT_MAX - 1  if n = 1
///        2) INT_MAX      if n is a prime
///        3) 0            if moebius(n) = 0
///        4) lpf - 1      if moebius(n) = 1
///        5) lpf          if moebius(n) = -1
///
///        factor[1] = (INT_MAX - 1) because 1 contributes to the
///        sum of the ordinary leaves S1(x, a) in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms.
///        The values above allow to replace the 1st if statement
///        below used in the S1(x, a) and S2(x, a) formulas by the
///        2nd new if statement which is obviously faster.
///
///        * Old: if (mu[n] != 0 && prime < lpf[n])
///        * New: if (prime < factor[n])
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdint.h>
#include <vector>

namespace primecount {

/// AbstractFactorTable contains static lookup tables
/// and is used to convert:
/// 1) A number into a FactorTable index
/// 2) A FactorTable index into a number
///
class AbstractFactorTable
{
public:
  static int64_t to_index(uint64_t number)
  {
    assert(number > 0);
    uint64_t q = number / 2310;
    uint64_t r = number % 2310;
    return 480 * q + coprime_indexes_[r];
  }

  static int64_t to_number(uint64_t index)
  {
    uint64_t q = index / 480;
    uint64_t r = index % 480;
    return 2310 * q + coprime_[r];
  }

  /// Returns the 1st number > 1 that is not divisible
  /// by 2, 3, 5, 7 and 11. Hence 13 is returned.
  ///
  static int64_t get_first_coprime()
  {
    return to_number(1);
  }

private:
  static const uint16_t coprime_[480];
  static const int16_t coprime_indexes_[2310];
};

template <typename T>
class FactorTable : public AbstractFactorTable
{
public:
  /// Factor numbers <= y
  FactorTable(int64_t y, int threads)
  {
    if (y > max())
      throw primecount_error("y must be <= FactorTable::max()");

    y = std::max<int64_t>(1, y);
    T T_MAX = std::numeric_limits<T>::max();
    factor_.resize(to_index(y) + 1, T_MAX);

    // mu(1) = 1.
    // 1 has zero prime factors, hence 1 has an even
    // number of prime factors. We use the least
    // significant bit to indicate whether the number
    // has an even or odd number of prime factors.
    factor_[0] ^= 1;

    int64_t sqrty = isqrt(y);
    int64_t thread_threshold = (int64_t) 1e7;
    threads = ideal_num_threads(threads, y, thread_threshold);
    int64_t thread_distance = ceil_div(y, threads);

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      int64_t low = 1;
      low += thread_distance * t;
      int64_t high = std::min(low + thread_distance, y);
      int64_t start = get_first_coprime() - 1;
      primesieve::iterator it(start);

      while (true)
      {
        int64_t i = 1;
        int64_t prime = it.next_prime();
        int64_t multiple = next_multiple(prime, low, &i);
        int64_t min_m = prime * get_first_coprime();

        if (min_m > high)
          break;

        for (; multiple <= high; multiple = prime * to_number(i++))
        {
          int64_t mi = to_index(multiple);
          // prime is smallest factor of multiple
          if (factor_[mi] == T_MAX)
            factor_[mi] = (T) prime;
          // the least significant bit indicates
          // whether multiple has an even (0) or odd (1)
          // number of prime factors
          else if (factor_[mi] != 0)
            factor_[mi] ^= 1;
        }

        if (prime <= sqrty)
        {
          int64_t j = 0;
          int64_t square = prime * prime;
          multiple = next_multiple(square, low, &j);

          // moebius(n) = 0
          for (; multiple <= high; multiple = square * to_number(j++))
            factor_[to_index(multiple)] = 0;
        }
      }
    }
  }

  /// mu_lpf(n) is a combination of the mu(n) (Möbius function)
  /// and lpf(n) (least prime factor) functions.
  /// mu_lpf(n) returns (with n = to_number(index)):
  ///
  /// 1) INT_MAX - 1  if n = 1
  /// 2) INT_MAX      if n is a prime
  /// 3) 0            if moebius(n) = 0
  /// 4) lpf - 1      if moebius(n) = 1
  /// 5) lpf          if moebius(n) = -1
  ///
  int64_t mu_lpf(int64_t index) const
  {
    return factor_[index];
  }

  /// Get the Möbius function value of the number
  /// n = to_number(index).
  ///
  /// https://en.wikipedia.org/wiki/Möbius_function
  /// mu(n) = 1 if n is a square-free integer with an even number of prime factors.
  /// mu(n) = −1 if n is a square-free integer with an odd number of prime factors.
  /// mu(n) = 0 if n has a squared prime factor.
  ///
  int64_t mu(int64_t index) const
  {
    if (factor_[index] == 0)
      return 0;
    else if (factor_[index] & 1)
      return -1;
    else
      return 1;
  }

  static maxint_t max()
  {
    maxint_t T_MAX = std::numeric_limits<T>::max();
    return ipow(T_MAX - 1, 2) - 1;
  }

private:

  /// Find the first multiple (of prime) > low which
  /// is not divisible by any prime <= 11.
  ///
  static int64_t next_multiple(int64_t prime,
                               int64_t low,
                               int64_t* index)
  {
    int64_t quotient = ceil_div(low, prime);
    int64_t i = std::max(*index, to_index(quotient));
    int64_t multiple = 0;

    for (; multiple <= low; i++)
      multiple = prime * to_number(i);

    *index = i;
    return multiple;
  }

  std::vector<T> factor_;
};

} // namespace

#endif
