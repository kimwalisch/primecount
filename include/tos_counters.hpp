///
/// @file   tos_counters.hpp
/// @brief  The counters data structure is a binary indexed
///         tree (a.k.a. Fenwick tree) that keeps track of
///         the number of unsieved elements (sieve[i] = 1)
///         in the sieve array. Whenever an element is
///         crossed-off for the first time in the sieve array
///         we update the counters data structure. Both
///         updating and querying the counters data structure
///         uses O(log n) operations.
///
///         The implementation is based on the paper:
///
///         Tomás Oliveira e Silva, Computing pi(x): the combinatorial method,
///         Revista do DETUA, vol. 4, no. 6, March 2006, pp. 767-768.
///         http://sweet.ua.pt/tos/bib/5.4.pdf
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef TOS_COUNTERS_HPP
#define TOS_COUNTERS_HPP

#include <stdint.h>
#include <vector>

namespace primecount {

/// Initialize the counters from the sieve array.
/// @pre segment_size is a power of 2
/// @pre sieve[i] = 1 for unsieved elements
///      sieve[i] = 0 for sieved elements
/// Runtime: O(N log N)
///
template <typename T1, typename T2>
inline void cnt_finit(const T1& sieve, std::vector<T2>& cnt, int64_t segment_size)
{
  segment_size >>= 1;
  cnt.resize(segment_size);

  for (int64_t i = 0; i < segment_size; i++)
  {
    cnt[i] = sieve[i * 2];
    int64_t k = (i + 1) & ~i;
    for (int64_t j = i; k >>= 1; j &= j - 1)
      cnt[i] += cnt[j - 1];
  }
}

/// Update (decrement by 1) the counters after
/// that an element has been crossed-off for the
/// first time in the sieve array.
/// Runtime: O(log N)
///
template <typename T>
inline void cnt_update(std::vector<T>& cnt, int64_t pos, int64_t segment_size)
{
  pos >>= 1;
  segment_size >>= 1;
  do
  {
    cnt[pos]--;
    pos |= pos + 1;
  }
  while (pos < segment_size);
}

/// Get the number of unsieved elements <= pos
/// in the current segment (sieve array).
/// Runtime: O(log N)
///
template <typename T>
inline int64_t get_sum(const std::vector<T>& cnt, int64_t pos)
{
  pos >>= 1;
  int64_t sum = cnt[pos++];
  for (; pos &= pos - 1; sum += cnt[pos - 1]);
  return sum;
}

} // namespace

#endif
