///
/// @file   moebius.cpp
/// @brief  Test the generate_moebius(n) function
/// @link   https://en.wikipedia.org/wiki/Moebius_function
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;
using namespace primecount;

vector<int> moebius =
{
    0, 1, -1, -1, 0, -1, 1, -1, 0, 0,
    1, -1, 0, -1, 1, 1, 0, -1, 0, -1,
    0, 1, 1, -1, 0, 0, 1, 0, 0, -1,
    -1, -1, 0, 1, 1, 1, 0, -1, 1, 1,
    0, -1, -1, -1, 0, 0, 1, -1, 0, 0,
    0, 1, 0, -1, 0, 1, 0, 1, 1, -1,
    0, -1, 1, 0, 0, 1, -1, -1, 0, 1,
    -1, -1, 0, -1, 1, 0, 0, 1, -1
};

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  auto mu = generate_moebius(moebius.size() - 1);

  for (size_t i = 1; i < mu.size(); i++)
  {
    cout << "mu(" << i << ") = " << mu[i];
    check(mu[i] == moebius[i]);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
