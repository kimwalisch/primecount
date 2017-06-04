///
/// @file   RiemannR.cpp
/// @brief  Test the Riemann R function
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

vector<int64_t> Ri_table =
{
                     4, // Ri(10^1)
                    25, // Ri(10^2)
                   168, // Ri(10^3)
                  1226, // Ri(10^4)
                  9587, // Ri(10^5)
                 78527, // Ri(10^6)
                664667, // Ri(10^7)
               5761551, // Ri(10^8)
              50847455, // Ri(10^9)
             455050683, // Ri(10^10)
          4118052494ll, // Ri(10^11)
         37607910542ll, // Ri(10^12)
        346065531065ll, // Ri(10^13)
       3204941731601ll, // Ri(10^14)
      29844570495886ll, // Ri(10^15)
     279238341360977ll, // Ri(10^16)
    2623557157055978ll  // Ri(10^17)
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
    size = Ri_table.size();

  for (size_t i = 0; i < size; i++)
  {
    int64_t x = ipow(10ll, i + 1);
    cout << "Ri(" << x << ") = " << Ri(x);
    check(Ri(x) == Ri_table[i]);
  }

  for (size_t i = 0; i < size; i++)
  {
    int64_t x = ipow(10ll, i + 1);
    cout << "Ri_inverse(" << Ri_table[i] << ") = " << Ri_inverse(Ri_table[i]);
    check(Ri_inverse(Ri_table[i]) < x &&
          Ri_inverse(Ri_table[i] + 1) >= x);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
