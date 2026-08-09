///
/// @file  FactorTableD.hpp
/// @brief The FactorTableD class combines the lpf[n] (least prime
///        factor), mpf[n] (max prime factor) and mu[n] (Möbius
///        function) lookup tables into a single factor[n] lookup
///        table which furthermore only contains entries for numbers
///        which are not divisible by 2, 3, 5, 7 and 11. The factor[n]
///        lookup table uses up to 28 times less memory than the
///        lpf[n], mpf[n] and mu[n] lookup tables! factor[n] uses only
///        2 bytes per entry for 32-bit numbers and 4 bytes per entry
///        for 64-bit numbers. The 16-bit ordinal format additionally
///        uses one Möbius bit per entry (about 2.125 bytes per entry).
///
///        The factor table concept was devised and implemented by
///        Christian Bau in 2003. Note that Tomás Oliveira e Silva
///        also suggests combining the mu[n] and lpf[n] lookup tables
///        in his paper. However Christian Bau's FactorTable data
///        structure uses only half as much memory and is also
///        slightly more efficient (uses fewer instructions) than the
///        data structure proposed by Tomás Oliveira e Silva.
///
///        By default we store the numerical least prime factor and
///        the Möbius sign in the factor[n] lookup table:
///
///        1) INT_MAX - 1  if n = 1
///        2) INT_MAX      if n is a prime
///        3) 0            if n has a prime factor > y
///        4) 0            if moebius(n) = 0
///        5) lpf - 1      if moebius(n) = 1
///        6) lpf          if moebius(n) = -1
///
///        factor[1] = (INT_MAX - 1) because 1 contributes to the
///        sum of the ordinary leaves S1(x, a) in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms.
///        The values above allow to replace the 1st if statement
///        below used in the D(x, y) formula by the 2nd new if
///        statement which is obviously faster.
///
///        * Old: if (mu[n] != 0 && lpf[n] > prime && mpf[n] <= y)
///        * New: if (factor[n] > prime)
///
///        FactorTableD<T, true> instead stores the ordinal of the
///        least prime factor. Its Möbius sign is stored in a separate
///        bitmap. This format allows a 16-bit factor table to be used
///        for much larger z values. In this mode the filter becomes:
///
///        * Ordinal: if (factor[n] > pi(prime))
///
///        In-depth description of the factor table data structure:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves-SIMD-Filtering.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLED_HPP
#define FACTORTABLED_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <BaseFactorTable.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <stdint.h>

namespace {

using namespace primecount;

template <typename T, bool store_ordinal = false>
class FactorTableD : public BaseFactorTable
{
public:
  using value_type = T;

  /// Factor numbers <= z
  FactorTableD(int64_t y,
               int64_t z,
               int threads)
  {
    if_unlikely(z > max())
      throw primecount_error("z must be <= FactorTable::max()");

    z = std::max<int64_t>(1, z);
    T T_MAX = pstd::numeric_limits<T>::max();
    factor_.resize(to_index(z) + 1);

    // mu(1) = 1.
    // 1 has zero prime factors, hence 1 has an even
    // number of prime factors. The numerical encoding uses the
    // least significant factor bit, the ordinal encoding uses its
    // separate Möbius bitmap.
    factor_[0] = T_MAX ^ 1;

    int64_t sqrtz = isqrt(z);
    int64_t thread_threshold = (int64_t) 1e7;
    threads = ideal_num_threads(z, threads, thread_threshold);
    int64_t thread_distance = ceil_div(z, threads);
    int64_t thread_alignment = coprime_indexes_.size();

    // In ordinal mode each thread must write to disjoint Möbius
    // bitmap words. Two wheel intervals contain 960 factor table
    // entries, which is divisible by 64.
    if (store_ordinal)
      thread_alignment *= 2;

    thread_distance += thread_alignment - thread_distance % thread_alignment;

    if (store_ordinal)
    {
      std::size_t words = (factor_.size() + 63) / 64;
      moebius_.resize(words);
      std::fill(moebius_.begin(), moebius_.end(), UINT64_MAX);
      moebius_[0] &= ~uint64_t(1);
    }

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      // Thread processes interval [low, high]
      int64_t low = thread_distance * t;
      int64_t high = low + thread_distance;
      low = std::max(first_coprime(), low + 1);
      high = std::min(high, z);

      if (low <= high)
      {
        // Default initialize memory to all bits set
        int64_t low_idx = to_index(low);
        int64_t size = (to_index(high) + 1) - low_idx;

        if (store_ordinal && t > 0)
          ASSERT(low_idx % 64 == 0);

        std::fill_n(&factor_[low_idx], size, T_MAX);

        int64_t start = first_coprime();
        int64_t stop = high / first_coprime();
        int64_t min_m = first_coprime() * first_coprime();
        primesieve::iterator it(start, stop);

        if (min_m <= high)
        {
          // pi(11) = 5, the first iterator prime is 13.
          int64_t prime_index = 5;

          while (true)
          {
            // Find multiples > prime
            int64_t i = 1;
            int64_t prime = it.next_prime();
            prime_index++;
            int64_t multiple = next_multiple(prime, low, &i);
            min_m = prime * first_coprime();

            if (min_m > high)
              break;

            for (; multiple <= high; multiple = prime * to_number(i++))
            {
              int64_t mi = to_index(multiple);
              // prime is the smallest factor of multiple
              if (factor_[mi] == T_MAX)
                factor_[mi] = encode_factor(prime, prime_index);
              // Toggle whether multiple has an even or odd number
              // of distinct prime factors.
              else if (factor_[mi] != 0)
                toggle_moebius(mi);
            }

            if (prime <= sqrtz)
            {
              int64_t j = 0;
              int64_t square = prime * prime;
              multiple = next_multiple(square, low, &j);

              // Sieve out numbers that are not square free
              // i.e. numbers for which moebius(n) = 0.
              for (; multiple <= high; multiple = square * to_number(j++))
                factor_[to_index(multiple)] = 0;
            }
          }
        }

        // Iterate over primes from [y+1, high]
        start = std::max(start, y + 1);

        if (start <= high)
        {
          it.jump_to(start, high);

          // y < prime <= z
          while (true)
          {
            int64_t i = 0;
            int64_t prime = it.next_prime();
            int64_t next = next_multiple(prime, low, &i);

            if (prime > high)
              break;

            // Sieve out primes > y &&
            // Sieve out numbers with prime factors > y
            for (; next <= high; next = prime * to_number(i++))
              factor_[to_index(next)] = 0;
          }
        }
      }
    }
  }

  /// Returns true if n (with n = to_number(index)) is a
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
  /// In ordinal mode composite leaves instead store pi(lpf),
  /// while the two sentinel values and zero keep the same meaning.
  ///
  int64_t is_leaf(int64_t index) const
  {
    return factor_[index];
  }

  const T* data() const
  {
    return factor_.data();
  }

  /// Convert the current prime into the value used by the hot
  /// factor_table[m] > filter comparison.
  int64_t get_filter_value(int64_t prime,
                           int64_t prime_index) const
  {
    return store_ordinal ? prime_index : prime;
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
    // mu(n) = 0 is disabled by default for performance
    // reasons, we only enable it for testing.
    #if defined(ENABLE_MU_0_TESTING)
      if (factor_[index] == 0)
        return 0;
    #else
      ASSERT(factor_[index] != 0);
    #endif

    if (store_ordinal)
    {
      std::size_t word = (std::size_t) index / 64;
      uint64_t bit = uint64_t(1) << (index % 64);
      return (moebius_[word] & bit) ? -1 : 1;
    }
    else if (factor_[index] & 1)
      return -1;
    else
      return 1;
  }

#if __cplusplus >= 201402L

  static constexpr int64_t max()
  {
    static_assert(sizeof(T) * 2 <= sizeof(uint64_t), "FactorTableD: sizeof(T) is too large!");
    static_assert(!store_ordinal ||
                  sizeof(T) <= sizeof(uint16_t),
                  "Ordinal FactorTableD only supports 8-bit and 16-bit values!");
    constexpr uint64_t MAX_T = pstd::numeric_limits<T>::max();
    constexpr uint64_t MAX_INT64_T = pstd::numeric_limits<int64_t>::max();
    constexpr uint64_t MAX_M_NUMERIC = (MAX_T - 1) * (MAX_T - 1) - 1;

    // Codes T_MAX - 1 and T_MAX are reserved for n = 1 and primes.
    // p(254) = 1609 and p(65534) = 821573 are therefore the first
    // unencodable least prime factors for 8-bit and 16-bit ordinals.
    constexpr uint64_t FIRST_UNENCODABLE_PRIME =
      sizeof(T) == sizeof(uint8_t) ? 1609 : 821573;
    constexpr uint64_t MAX_M_ORDINAL =
      FIRST_UNENCODABLE_PRIME * FIRST_UNENCODABLE_PRIME - 1;
    constexpr uint64_t MAX_M = store_ordinal ? MAX_M_ORDINAL : MAX_M_NUMERIC;
    return (int64_t) std::min(MAX_M, MAX_INT64_T);
  }

#else

  static int64_t max()
  {
    static_assert(sizeof(T) * 2 <= sizeof(uint64_t), "FactorTableD: sizeof(T) is too large!");
    static_assert(!store_ordinal ||
                  sizeof(T) <= sizeof(uint16_t),
                  "Ordinal FactorTableD only supports 8-bit and 16-bit values!");
    uint64_t MAX_T = pstd::numeric_limits<T>::max();
    uint64_t MAX_INT64_T = pstd::numeric_limits<int64_t>::max();
    uint64_t MAX_M_NUMERIC = (MAX_T - 1) * (MAX_T - 1) - 1;
    uint64_t FIRST_UNENCODABLE_PRIME =
      sizeof(T) == sizeof(uint8_t) ? 1609 : 821573;
    uint64_t MAX_M_ORDINAL =
      FIRST_UNENCODABLE_PRIME * FIRST_UNENCODABLE_PRIME - 1;
    uint64_t MAX_M = store_ordinal ? MAX_M_ORDINAL : MAX_M_NUMERIC;
    return (int64_t) std::min(MAX_M, MAX_INT64_T);
  }

#endif

private:
  static T encode_factor(int64_t prime,
                         int64_t prime_index)
  {
    if (store_ordinal)
    {
      ASSERT(prime_index <= (int64_t) pstd::numeric_limits<T>::max() - 2);
      return (T) prime_index;
    }
    else
      return (T) prime;
  }

  void toggle_moebius(int64_t index)
  {
    if (store_ordinal)
    {
      std::size_t word = (std::size_t) index / 64;
      uint64_t bit = uint64_t(1) << (index % 64);
      moebius_[word] ^= bit;
    }
    else
      factor_[index] ^= 1;
  }

  Vector<T> factor_;
  Vector<uint64_t> moebius_;
};

template <typename T>
using FactorTableDOrdinal = FactorTableD<T, true>;

} // namespace

#endif
