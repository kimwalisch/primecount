///
/// @file  SemiprimeEventSorter.hpp
/// @brief Sort distinct semiprime events by decreasing n using two
///        wheel-compressed bitmaps for presence and weight.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SEMIPRIME_EVENT_SORTER_HPP
#define SEMIPRIME_EVENT_SORTER_HPP

#include <BitSieve240.hpp>
#include <ctz.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

class SemiprimeEventSorter : private BitSieve240
{
public:
  /// Events must have distinct n=q*r in (lower, upper], where p<=q<=r
  /// are primes, and weight 1 or 2. Retain one sorter per D work unit
  /// so its scratch allocations are reused between segments and primes.
  template <typename Event>
  void sort(Vector<Event>& events, int64_t p, int64_t lower, int64_t upper)
  {
    if (events.size() < 2)
      return;

    ASSERT(0 <= lower && lower < upper);
    uint64_t first = uint64_t(lower + 1) / 240;
    uint64_t last = uint64_t(upper) / 240;
    uint64_t words64 = last - first + 1;
    std::size_t count = events.size();

    // Check the full interval before narrowing on 32-bit platforms.
    // Compare bitmap words with N * floor(log2(N)) without multiplying.
    // Small primes may produce events that are not coprime to 30.
    if (p <= 5 ||
        words64 > SIZE_MAX ||
        (words64 - 1) / count >= ilog2(count))
    {
      std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        return a.n > b.n;
      });
      return;
    }

    std::size_t words = std::size_t(words64);
    present_.resize(words);
    weight1_.resize(words);
    std::fill(present_.begin(), present_.end(), 0);
    std::fill(weight1_.begin(), weight1_.end(), 0);

    for (const Event& event : events)
    {
      ASSERT(lower < event.n && event.n <= upper);
      ASSERT(event.weight == 1 || event.weight == 2);
      uint64_t n = uint64_t(event.n);
      uint64_t word = n / 240 - first;
      uint64_t residue = n % 240;
      ASSERT(residue != 0);

      // Reverse the bits so ctz64 enumerates n in decreasing order
      uint64_t bit = set_bit_[240 - residue];
      ASSERT(bit != 0);
      ASSERT((present_[word] & bit) == 0);
      present_[word] |= bit;
      if (event.weight == 1)
        weight1_[word] |= bit;
    }

    std::size_t written = 0;
    for (std::size_t word = words; word-- > 0;)
    {
      uint64_t bits = present_[word];
      uint64_t weight1 = weight1_[word];
      uint64_t base = (first + word) * 240;
      while (bits)
      {
        int bit = ctz64(bits);
        int64_t n = int64_t(base + bit_values_[63 - bit]);
        int64_t weight = 2 - ((weight1 >> bit) & 1);
        events[written++] = Event{n, weight};
        bits &= bits - 1;
      }
    }
    ASSERT(written == events.size());
  }

private:
  Vector<uint64_t> present_;
  Vector<uint64_t> weight1_;
};

} // namespace

#endif
