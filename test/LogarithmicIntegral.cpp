///
/// @file   LogarithmicIntegral.cpp
/// @brief  Test the offset logarithmic integral function
///         Li(x) = li(x) - li(2)
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;
using namespace primecount;

vector<int64_t> Li_table =
{
                     5, // Li(10^1)
                    29, // Li(10^2)
                   176, // Li(10^3)
                  1245, // Li(10^4)
                  9628, // Li(10^5)
                 78626, // Li(10^6)
                664917, // Li(10^7)
               5762208, // Li(10^8)
              50849233, // Li(10^9)
             455055613, // Li(10^10)
          4118066399ll, // Li(10^11)
         37607950279ll, // Li(10^12)
        346065645809ll, // Li(10^13)
       3204942065690ll, // Li(10^14)
      29844571475286ll, // Li(10^15)
     279238344248555ll, // Li(10^16)
    2623557165610820ll  // Li(10^17)
};

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  size_t size = 15;

  if (sizeof(long double) > 8)
    size = Li_table.size();

  for (size_t i = 0; i < size; i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    cout << "Li(" << x << ") = " << Li(x);
    check(Li(x) == Li_table[i]);
  }

  for (size_t i = 0; i < size; i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    cout << "Li_inverse(" << Li_table[i] << ") = " << Li_inverse(Li_table[i]);
    check(Li_inverse(Li_table[i]) <= x &&
          Li_inverse(Li_table[i] + 1) > x);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
