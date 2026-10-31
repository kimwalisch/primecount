///
/// @file  SemiprimeEventSorter.hpp
/// @brief Sort the distinct composite events of hyperbola_sum() by
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
#include <popcnt.hpp>
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
  /// Events must have distinct n in (min_n, max_n],
  /// coprime to 30, and weight 1 or 2. Retain
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
    uint64_t word_count = last - first + 1;
    uint64_t count = events.size();

    // Fall back to std::sort if the bitmap has more than
    // N * floor(log2(N)) words, which is checked without
    // multiplying.
    if ((word_count - 1) / count >= ilog2(count))
    {
      std::sort(events.begin(), events.end(), [](const SemiprimeEvent& a,
                                                 const SemiprimeEvent& b) {
        return a.n > b.n;
      });
      return;
    }

    present_.clear();
    weight1_.clear();

    present_.resize(word_count);
    weight1_.resize(word_count);

    std::fill(present_.begin(), present_.end(), 0);
    std::fill(weight1_.begin(), weight1_.end(), 0);

    // Insert the events into two bitmaps: present_ marks all
    // events and weight1_ marks the events with weight 1.
    // Since all events are coprime to 30, each 64-bit word
    // covers 240 integers. The bitmaps store offset = high - n
    // instead of n, so that increasing offset (word and bit)
    // corresponds to decreasing n.
    uint64_t high = (last + 1) * 240;

    for (const SemiprimeEvent& event : events)
    {
      uint64_t n = event.n;
      ASSERT(min_n < n && n <= max_n);
      ASSERT(event.weight == 1 || event.weight == 2);
      uint64_t offset = high - n;
      uint64_t word = offset / 240;
      uint64_t bit = set_bit_[offset % 240];
      ASSERT(bit != 0);
      ASSERT((present_[word] & bit) == 0);
      present_[word] |= bit;
      if (event.weight == 1)
        weight1_[word] |= bit;
    }

    // Extra space for the unconditional writes below
    events.resize(count + 4);
    uint64_t written = 0;

    // Read the bitmaps word by word in increasing order of
    // offset, i.e. in decreasing order of n, and write the
    // events back into the events vector in sorted order.
    for (uint64_t word = 0; word < word_count; word++)
    {
      uint64_t bits = present_[word];
      uint64_t bit_count = popcnt64(bits);
      uint64_t weight1 = weight1_[word];
      uint64_t base = high - word * 240;
      SemiprimeEvent* out = &events[written];
      written += bit_count;

      // Decode the events of the current word 4 at a time,
      // from the lowest to the highest set bit. The last
      // iteration may decode up to 3 dummy events, which are
      // later overwritten or discarded. Decoding 4 events
      // unconditionally avoids most branch mispredictions.
      NO_UNROLL_LOOP
      NO_VECTORIZE_LOOP
      for (uint64_t i = 0; i < bit_count; i += 4)
      {
        out[i + 0] = next_event(bits, base, weight1); bits &= bits - 1;
        out[i + 1] = next_event(bits, base, weight1); bits &= bits - 1;
        out[i + 2] = next_event(bits, base, weight1); bits &= bits - 1;
        out[i + 3] = next_event(bits, base, weight1); bits &= bits - 1;
      }
    }

    events.resize(count);
    ASSERT(written == events.size());
  }

private:
  /// Decode the event of the lowest set bit.
  /// If bits = 0, a dummy event is returned.
  ALWAYS_INLINE static SemiprimeEvent
  next_event(uint64_t bits,
             uint64_t base,
             uint64_t weight1)
  {
    int bit_index = ctz64(bits | (1ull << 63));
    int64_t n = base - bit_values_[bit_index];
    int64_t weight = 2 - ((weight1 >> bit_index) & 1);
    return SemiprimeEvent{n, weight};
  }

  Vector<uint64_t> present_;
  Vector<uint64_t> weight1_;
};

} // namespace

#endif
