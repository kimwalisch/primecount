///
/// @file   Ri.cpp
/// @brief  Test the Riemann R function.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>

using std::log;
using std::max;
using std::size_t;
using namespace primecount;

std::vector<int64_t> Ri_table =
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
  3204941731601ll  // Ri(10^14)
};

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  for (size_t i = 0; i < Ri_table.size(); i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    std::cout << "Ri(" << x << ") = " << Ri(x);
    check(Ri(x) == Ri_table[i]);
  }

  for (size_t i = 0; i < Ri_table.size(); i++)
  {
    int p = (int) i + 1;
    int64_t x = ipow(10ll, p);
    std::cout << "Ri_inverse(" << Ri_table[i] << ") = " << Ri_inverse(Ri_table[i]);
    check(Ri_inverse(Ri_table[i]) < x &&
          Ri_inverse(Ri_table[i] + 1) >= x);
  }

  int64_t x;

  // Sanity checks for tiny values of Ri(x)
  for (x = 0; x < 10000; x++)
  {
    int64_t rix = Ri(x);
    double logx = log(max((double) x, 2.0));

    if (rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "Ri(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of Ri(x)
  for (; x < 100000; x += 101)
  {
    int64_t rix = Ri(x);
    double logx = log(max((double) x, 2.0));

    if (rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "Ri(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for tiny values of Ri_inverse(x)
  for (x = 2; x < 1000; x++)
  {
    int64_t res = Ri_inverse(x);
    double logx = log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "Ri_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of Ri_inverse(x)
  for (; x < 100000; x += 101)
  {
    int64_t res = Ri_inverse(x);
    double logx = log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "Ri_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
