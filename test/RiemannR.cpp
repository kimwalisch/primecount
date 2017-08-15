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
  29844570495886ll  // Ri(10^15)
};

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  for (size_t i = 0; i < Ri_table.size(); i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    cout << "Ri(" << x << ") = " << Ri(x);
    check(Ri(x) == Ri_table[i]);
  }

  for (size_t i = 0; i < Ri_table.size(); i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    cout << "Ri_inverse(" << Ri_table[i] << ") = " << Ri_inverse(Ri_table[i]);
    check(Ri_inverse(Ri_table[i]) < x &&
          Ri_inverse(Ri_table[i] + 1) >= x);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
