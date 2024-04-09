///
/// @file   Li.cpp
/// @brief  Test the Eulerian logarithmic integral function.
///         Li(x) = li(x) - li(2)
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <calculator.hpp>
#include <int128_t.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using std::size_t;
using namespace primecount;

/// Generated using Mathematica:
/// Table[IntegerPart[li[k]-li[2]], {k, 2, 99}]
std::array<int64_t, 100> Li_tiny =
{
  0, 0, 0, 1, 1, 2, 3, 3, 4, 4,
  5, 5, 5, 6, 6, 7, 7, 7, 8, 8,
  8, 9, 9, 9, 10, 10, 10, 11, 11, 11,
  11, 12, 12, 12, 13, 13, 13, 13, 14, 14,
  14, 15, 15, 15, 15, 16, 16, 16, 16, 17,
  17, 17, 17, 18, 18, 18, 18, 19, 19, 19,
  19, 20, 20, 20, 20, 21, 21, 21, 21, 22,
  22, 22, 22, 23, 23, 23, 23, 23, 24, 24,
  24, 24, 25, 25, 25, 25, 25, 26, 26, 26,
  26, 27, 27, 27, 27, 27, 28, 28, 28, 28
};

std::vector<int64_t> Li_table =
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
 3204942065690ll  // Li(10^14)
};

#if defined(HAVE_FLOAT128)

std::vector<std::string> Li_f128 =
{
                 "29844571475286", // Li(10^15)
                "279238344248555", // Li(10^16)
               "2623557165610820", // Li(10^17)
              "24739954309690413", // Li(10^18)
             "234057667376222381", // Li(10^19)
            "2220819602783663482", // Li(20^20)
           "21127269486616126181", // Li(20^21)
          "201467286691248261497", // Li(20^22)
         "1925320391614054155137", // Li(20^23)
        "18435599767366347775143", // Li(20^24)
       "176846309399198930392618", // Li(20^25)
      "1699246750872593033005722", // Li(20^26)
     "16352460426842189113085404", // Li(20^27)
    "157589269275974838158399970", // Li(20^28)
   "1520698109714276717287880526", // Li(20^29)
  "14692398897720447639079087668"  // Li(20^30)
};

#endif

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  for (int64_t x = 0; x < (int64_t) Li_tiny.size(); x++)
  {
    std::cout << "Li(" << x << ") = " << Li(x);
    check(Li(x) == Li_tiny[x]);
  }

  {
    int64_t x = 10;
    for (size_t i = 0; i < Li_table.size(); i++)
    {
      std::cout << "Li(" << x << ") = " << Li(x);
      check(Li(x) == Li_table[i]);
      x *= 10;
    }
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  {
    int128_t x = ipow<15>((int128_t) 10);
    for (size_t i = 0; i < Li_f128.size(); i++)
    {
      std::ostringstream oss;
      oss << Li(x);
      std::cout << "Li(" << x << ") = " << oss.str();
      check(oss.str() == Li_f128[i]);
      x *= 10;
    }
  }

#endif

  for (int64_t x = 1; x < (int64_t) Li_tiny.size(); x++)
  {
    int64_t y = Li_tiny[x];
    std::cout << "Li_inverse(" << y << ") = " << Li_inverse(y);
    check(Li_inverse(y) < x &&
          Li_inverse(y + 1) >= x);
  }

  {
    int64_t x = 10;
    for (size_t i = 0; i < Li_table.size(); i++)
    {
      int64_t y = Li_table[i];
      std::cout << "Li_inverse(" << y << ") = " << Li_inverse(y);
      check(Li_inverse(y) < x &&
            Li_inverse(y + 1) >= x);
      x *= 10;
    }
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  {
    int128_t x = ipow<15>((int128_t) 10);
    for (size_t i = 0; i < Li_f128.size(); i++)
    {
      int128_t y = calculator::eval<int128_t>(Li_f128[i]);
      std::cout << "Li_inverse(" << y << ") = " << Li_inverse(y);
      check(Li_inverse(y) < x &&
            Li_inverse(y + 1) >= x);
      x *= 10;
    }
  }

#endif

  {
    // Li(9760) = 1219.000098
    int64_t x = 9760;
    int64_t y = 1219;
    std::cout << "Li(" << x << ") = " << Li(x);
    check(Li(x) == y);

    std::cout << "Li_inverse(" << y << ") = " << Li_inverse(y);
    check(Li_inverse(y) < x &&
          Li_inverse(y + 1) >= x);
  }

  {
    // Li(9494) = 1189.9997
    int64_t x = 9494;
    int64_t y = 1189;
    std::cout << "Li(" << x << ") = " << Li(x);
    check(Li(x) == y);

    std::cout << "Li_inverse(" << y << ") = " << Li_inverse(y);
    check(Li_inverse(y) < x &&
          Li_inverse(y + 1) >= x);
  }

  // Sanity checks for tiny values of Li(x)
  int64_t x;
  for (x = 0; x < 10000; x++)
  {
    int64_t Lix = Li(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (Lix < 0 ||
        (x >= 11 && Lix < x / logx) ||
        (x >= 2  && Lix > x * logx))
    {
      std::cout << "Li(" << x << ") = " << Lix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of Li(x)
  for (; x < 100000; x += 101)
  {
    int64_t Lix = Li(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (Lix < 0 ||
        (x >= 11 && Lix < x / logx) ||
        (x >= 2  && Lix > x * logx))
    {
      std::cout << "Li(" << x << ") = " << Lix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for tiny values of Li_inverse(x)
  for (x = 2; x < 1000; x++)
  {
    int64_t res = Li_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 4 && res > x * logx * logx))
    {
      std::cout << "Li_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of Li_inverse(x)
  for (; x < 100000; x += 101)
  {
    int64_t res = Li_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 4 && res > x * logx * logx))
    {
      std::cout << "Li_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  {
    int64_t x = port::numeric_limits<int64_t>::max() / 10;
    int64_t res = Li_inverse(x);
    if (res != port::numeric_limits<int64_t>::max())
    {
      std::cout << "Li_inverse(" << x << ") != INT64_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }

#if defined(HAVE_INT128_T)
  {
    int128_t x = port::numeric_limits<int128_t>::max() / 10;
    int128_t res = Li_inverse(x);
    if (res != port::numeric_limits<int128_t>::max())
    {
      std::cout << "Li_inverse(" << x << ") != INT128_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
