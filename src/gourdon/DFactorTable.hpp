///
/// @file  DFactorTable.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef DFACTORTABLE_HPP
#define DFACTORTABLE_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <FactorTable.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdint.h>
#include <vector>

namespace primecount {

template <typename T>
class DFactorTable : public AbstractFactorTable
{
public:
  /// Factor numbers <= z
  DFactorTable(int64_t y, 
               int64_t z, 
               int threads)
  {
    if (z > max())
      throw primecount_error("z must be <= FactorTable::max()");

    z = std::max<int64_t>(1, z);
    T T_MAX = std::numeric_limits<T>::max();
    factor_.resize(get_index(z) + 1, T_MAX);

    // mu(1) = 1.
    // 1 has zero prime factors, hence 1 has an even
    // number of prime factors. We use the least
    // significant bit to indicate whether the number
    // has an even or odd number of prime factors.
    factor_[0] ^= 1;

    int64_t sqrtz = isqrt(z);
    int64_t thread_threshold = ipow(10, 7);
    threads = ideal_num_threads(threads, z, thread_threshold);
    int64_t thread_distance = ceil_div(z, threads);

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      int64_t low = 1;
      low += thread_distance * t;
      int64_t high = std::min(low + thread_distance, z);
      int64_t start = get_first_coprime() - 1;
      primesieve::iterator it(start);

      while (true)
      {
        // Find multiples > prime
        int64_t i = 1;
        int64_t prime = it.next_prime();
        int64_t multiple = next_multiple(prime, low, &i);
        int64_t min_m = prime * get_first_coprime();

        if (min_m > high)
          break;

        // Find the smallest prime factors of the
        // integers inside ]low, high].
        for (; multiple <= high; multiple = prime * get_number(i++))
        {
          int64_t mi = get_index(multiple);
          // prime is the smallest factor of multiple
          if (factor_[mi] == T_MAX)
            factor_[mi] = (T) prime;
          // the least significant bit indicates
          // whether multiple has an even (0) or odd (1)
          // number of prime factors
          else if (factor_[mi] != 0)
            factor_[mi] ^= 1;
        }

        if (prime <= sqrtz)
        {
          int64_t j = 0;
          int64_t square = prime * prime;
          multiple = next_multiple(square, low, &j);

          // Sieve out numbers that are not square free
          // i.e. numbers for which moebius(n) = 0.
          for (; multiple <= high; multiple = square * get_number(j++))
            factor_[get_index(multiple)] = 0;
        }
      }

      start = std::max(start, y);
      it.skipto(start, z);

      while (true)
      {
        int64_t i = 0;
        int64_t prime = it.next_prime();
        int64_t next = next_multiple(prime, low, &i);

        if (prime > high)
          break;

        // Sieve out primes > y &&
        // Sieve out numbers with prime factors > y
        for (; next <= high; next = prime * get_number(i++))
          factor_[get_index(next)] = 0;
      }
    }
  }

  /// Returns true if n (with n = get_number(index)) is a
  /// hard special leaf in the D formula of Xavier
  /// Gourdon's prime counting algorithm.
  ///
  /// Return value:
  ///
  /// 1) INT_MAX - 1  if n = 1
  /// 2) INT_MAX      if n is a prime
  /// 3) 0            if n has a prime factor > y
  /// 4) 0            if moebius(n) = 0
  /// 5) lpf - 1      if moebius(n) = 1
  /// 6) lpf          if moebius(n) = -1
  ///
  int64_t is_leaf(int64_t index) const
  {
    return factor_[index];
  }

  /// Get the Möbius function value of the number
  /// n = get_number(index).
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
    int64_t i = std::max(*index, get_index(quotient));
    int64_t multiple = 0;

    for (; multiple <= low; i++)
      multiple = prime * get_number(i);

    *index = i;
    return multiple;
  }

  std::vector<T> factor_;
};

} // namespace

#endif
