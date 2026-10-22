///
/// @file  SemiprimeEventSorter.cpp
/// @brief Test semiprime event ordering and weight preservation.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <SemiprimeEventSorter.hpp>
#include <Vector.hpp>
#include <primesieve/iterator.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <cstdlib>

using namespace primecount;

struct Event
{
  int64_t n;
  int64_t weight;
};

bool is_prime(int64_t n)
{
  for (int64_t d = 2; d <= n / d; d++)
    if (n % d == 0)
      return false;
  return n >= 2;
}

Vector<Event> trial_events(int64_t p, int64_t lower, int64_t upper)
{
  Vector<Event> events;
  for (int64_t n = lower + 1; n <= upper; n++)
    for (int64_t q = 2; q <= n / q; q++)
      if (n % q == 0)
      {
        int64_t r = n / q;
        if (q >= p && is_prime(r))
          events.push_back({n, q == p || q == r ? 1 : 2});
        break;
      }
  return events;
}

void check(SemiprimeEventBitmap& sorter, Vector<Event>& events,
           int64_t p, int64_t lower, int64_t upper)
{
  Vector<Event> expected(events.size());
  if (!events.empty())
  {
    std::copy(events.begin(), events.end(), expected.begin());
    std::sort(expected.begin(), expected.end(), [](const Event& a, const Event& b) {
      return a.n > b.n;
    });
    std::reverse(events.begin(), events.end());
    std::rotate(events.begin(), events.begin() + events.size() / 3, events.end());
  }
  sorter.sort(events, p, lower, upper);

  if (events.size() != expected.size())
  {
    std::cerr << "Semiprime event count changed\n";
    std::exit(1);
  }

  for (std::size_t i = 0; i < expected.size(); i++)
    if (events[i].n != expected[i].n || events[i].weight != expected[i].weight)
    {
      std::cerr << "Semiprime events differ: p=" << p
                << ", lower=" << lower << ", upper=" << upper
                << ", index=" << i << '\n';
      std::exit(1);
    }
}

int main()
{
  SemiprimeEventBitmap sorter;

  // Trial division independently supplies squares and both event weights
  for (int64_t p : {2, 3, 5, 7, 23})
  {
    auto events = trial_events(p, 0, 5000);
    check(sorter, events, p, 0, 5000);
  }

  const Array<int64_t, 9> endpoints = {0, 48, 49, 239, 240, 241, 666, 840, 898};
  for (int64_t lower : endpoints)
    for (int64_t width : {0, 1, 2, 239, 240, 241})
    {
      auto events = trial_events(7, lower, lower + width);
      check(sorter, events, 7, lower, lower + width);
    }

  // Exercise high bitmap addresses without allocating a large prime table
  const Array<int64_t, 2> starts = {4294963200, 1000000000000};
  const Array<int64_t, 8> factors = {7, 11, 13, 17, 19, 23, 29, 31};
  for (int64_t lower : starts)
  {
    int64_t upper = lower + 8192;
    Vector<Event> events;
    for (int64_t q : factors)
    {
      primesieve::iterator iter(lower / q + 1, upper / q);
      for (int64_t r = iter.next_prime(); r <= upper / q; r = iter.next_prime())
        events.push_back({q * r, q == 7 ? 1 : 2});
    }
    check(sorter, events, 7, lower, upper);

    // The same sparse events over a huge interval must avoid a huge bitmap
    check(sorter, events, 7, 0, upper);
    check(sorter, events, 7, 0, 240 * (1ull << 32));
  }

  // Reuse existing allocations with changing partial first and last blocks
  for (int64_t lower = 529; lower < 5000; lower += 137)
  {
    auto events = trial_events(23, lower, lower + 997);
    check(sorter, events, 23, lower, lower + 997);
  }

  std::cout << "All tests passed successfully!" << std::endl;
  return 0;
}
