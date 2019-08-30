///
/// @file   segmented_pi_table.cpp
/// @brief  Test the SegmentedPiTable class
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount.hpp>
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
  int threads = get_num_threads();

  PiTable pi(limit);
  SegmentedPiTable segmentedPi(limit, segment_size, threads);
  int64_t i = 0;

  // Check small pi(x) values
  for (; i <= 1000; i++)
  {
    while (i >= segmentedPi.high())
      segmentedPi.next();

    cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  // Check large pi(x) values
  for (; i < limit; i += dist2(gen))
  {
    while (i >= segmentedPi.high())
      segmentedPi.next();

    cout << "segmentedPi(" << i << ") = " << segmentedPi[i];
    check(segmentedPi[i] == pi[i]);
  }

  while (limit >= segmentedPi.high())
    segmentedPi.next();

  // Check max pi(x) value
  cout << "segmentedPi(" << limit << ") = " << segmentedPi[limit];
  check(segmentedPi[limit] == pi[limit]);

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
