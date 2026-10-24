///
/// @file  SemiprimeEventSorter.hpp
/// @brief Sort the distinct semiprime events of hyperbola_sum() by
///        decreasing n using two wheel-compressed bitmaps for presence
///        and weight, or std::sort if the bitmaps would be too large.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
/// Copyright (C) 2026 Chopper Lin, <choppers@126.com>
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

struct SemiprimeEvent
{
  int64_t n;
  int64_t weight;
};

class SemiprimeEventSorter : private BitSieve240
{
public:
  /// Events must have distinct n = q*r in (min_n, max_n],
  /// where 5 < q <= r are primes, and weight 1 or 2. Retain
  /// one sorter per D work unit so its scratch allocations
  /// are reused between segments and primes.
  void sort(Vector<SemiprimeEvent>& events,
            uint64_t min_n,
            uint64_t max_n)
  {
    if (events.size() < 2)
      return;

    ASSERT(min_n < max_n);
    uint64_t first = (min_n + 1) / 240;
    uint64_t last = max_n / 240;
    uint64_t words64 = last - first + 1;
    uint64_t count = events.size();

    // Fall back to std::sort if the bitmap has more than
    // N * floor(log2(N)) words, which is checked without
    // multiplying.
    if ((words64 - 1) / count >= ilog2(count))
    {
      std::sort(events.begin(), events.end(), [](const SemiprimeEvent& a,
                                                 const SemiprimeEvent& b) {
        return a.n > b.n;
      });
      return;
    }

    present_.clear();
    weight1_.clear();
    present_.resize(words64);
    weight1_.resize(words64);

    std::fill(present_.begin(), present_.end(), 0);
    std::fill(weight1_.begin(), weight1_.end(), 0);

    for (const SemiprimeEvent& event : events)
    {
      uint64_t n = event.n;
      ASSERT(min_n < n && n <= max_n);
      ASSERT(event.weight == 1 || event.weight == 2);
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

    uint64_t written = 0;

    for (uint64_t word = words64; word-- > 0;)
    {
      uint64_t bits = present_[word];
      uint64_t weight1 = weight1_[word];
      uint64_t base = (first + word) * 240;

      while (bits)
      {
        int bit = ctz64(bits);
        int64_t n = base + bit_values_[63 - bit];
        int64_t weight = 2 - ((weight1 >> bit) & 1);
        events[written++] = SemiprimeEvent{n, weight};
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
