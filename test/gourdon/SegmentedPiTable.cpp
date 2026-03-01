///
/// @file   SegmentedPiTable.cpp
/// @brief  Test the SegmentedPiTable class
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <random>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(9000000, 10000000);
  std::uniform_int_distribution<int> dist2(1, 1000);

  int64_t limit = dist(gen);
  int64_t segment_size = isqrt(limit);
  segment_size += 128 - segment_size % 128;
  int threads = 1;

  PiTable pi(limit, threads);
  SegmentedPiTable segmentedPi;

  int64_t i = 0;
  int64_t low = 0;
  int64_t high = segment_size;
  segmentedPi.init(low, high);

  // Check small pi(x) values
  for (; i <= 1000; i++)
  {
    while (high <= i)
    {
      low = high;
      high = low + segment_size;
      segmentedPi.init(low, high);
    }

    std::cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  // Check large pi(x) values
  for (; i < limit; i += dist2(gen))
  {
    while (high <= i)
    {
      low = high;
      high = low + segment_size;
      segmentedPi.init(low, high);
    }

    std::cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  while (high < limit)
  {
    low = high;
    high = low + segment_size;
    segmentedPi.init(low, high);
  }

  // Check max pi(x) value.
  // PiTable can lookup numbers <= limit.
  // SegmentedPiTable can lookup numbers < limit.
  std::cout << "segmentedPi(" << limit-1 << ") = " << segmentedPi[limit-1];
  check(segmentedPi[limit-1] == pi[limit-1]);

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
