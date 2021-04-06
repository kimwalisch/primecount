///
/// @file   segmented_pi_table.cpp
/// @brief  Test the SegmentedPiTable class
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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

using namespace std;
using namespace primecount;

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(90000000, 100000000);
  uniform_int_distribution<int> dist2(1, 1000);

  int64_t limit = dist(gen);
  int64_t segment_size = iroot<3>(limit);
  segment_size += 240 - segment_size % 240;
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
    while (i >= high)
    {
      low = high;
      high = low + segment_size;
      segmentedPi.init(low, high);
    }

    cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  // Check large pi(x) values
  for (; i < limit; i += dist2(gen))
  {
    while (i >= high)
    {
      low = high;
      high = low + segment_size;
      segmentedPi.init(low, high);
    }

    cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  while (limit > high)
  {
    low = high;
    high = low + segment_size;
    segmentedPi.init(low, high);
  }

  // Check max pi(x) value.
  // PiTable can lookup numbers <= limit.
  // SegmentedPiTable can lookup numbers < limit.
  cout << "segmentedPi(" << limit-1 << ") = " << segmentedPi[limit-1];
  check(segmentedPi[limit-1] == pi[limit-1]);

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
