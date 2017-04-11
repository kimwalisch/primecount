///
/// @file  libdivide.cpp
/// @brief Test the branchfree divider of libdivide.h
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <libdivide.h>
#include <stdint.h>
#include <limits>

uint64_t dividends[20] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 63, 101, 511,
    1 << 5, 1 << 9, 1 << 20,
    std::numeric_limits<int32_t>::max(),
    std::numeric_limits<uint32_t>::max(),
    std::numeric_limits<int64_t>::max(),
    std::numeric_limits<uint64_t>::max()
};

int main()
{
  for (int i = 0; i < 20; i++)
  {
    for (int j = 2; j < 10000; j++)
    {
      libdivide::divider<uint64_t, libdivide::BRANCHFREE> divider(j);
      if (dividends[i] / j != dividends[i] / divider)
        return 1;
    }
  }

  for (int i = 0; i < 20; i++)
  {
    for (int j = 2; j < 20; j++)
    {
      libdivide::divider<uint64_t, libdivide::BRANCHFREE> divider(dividends[j]);
      if (dividends[i] / dividends[j] != dividends[i] / divider)
        return 1;
    }
  }

  return 0;
}
