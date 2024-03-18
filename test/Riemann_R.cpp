///
/// @file   Riemann_R.cpp
/// @brief  Test the Riemann R function.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <calculator.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using std::size_t;
using namespace primecount;

/// Generated using Mathematica:
/// Table[IntegerPart[RiemannR[k]], {k, 0, 99}]
std::array<int64_t, 100> RiemannR_tiny =
{
  0, 1, 1, 2, 2, 2, 3, 3, 3, 4,
  4, 4, 5, 5, 5, 6, 6, 6, 6, 7,
  7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
  10, 10, 10, 10, 11, 11, 11, 11, 12, 12,
  12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
  14, 15, 15, 15, 15, 16, 16, 16, 16, 17,
  17, 17, 17, 17, 18, 18, 18, 18, 18, 19,
  19, 19, 19, 20, 20, 20, 20, 20, 21, 21,
  21, 21, 21, 22, 22, 22, 22, 23, 23, 23,
  23, 23, 24, 24, 24, 24, 24, 25, 25, 25
};

std::vector<int64_t> RiemannR_table =
{
               4, // RiemannR(10^1)
              25, // RiemannR(10^2)
             168, // RiemannR(10^3)
            1226, // RiemannR(10^4)
            9587, // RiemannR(10^5)
           78527, // RiemannR(10^6)
          664667, // RiemannR(10^7)
         5761551, // RiemannR(10^8)
        50847455, // RiemannR(10^9)
       455050683, // RiemannR(10^10)
    4118052494ll, // RiemannR(10^11)
   37607910542ll, // RiemannR(10^12)
  346065531065ll, // RiemannR(10^13)
 3204941731601ll  // RiemannR(10^14)
};

#if defined(HAVE_FLOAT128)

std::vector<std::string> RiemannR_f128 =
{
                 "29844570495886", // RiemannR(10^15)
                "279238341360977", // RiemannR(10^16)
               "2623557157055978", // RiemannR(10^17)
              "24739954284239494", // RiemannR(10^18)
             "234057667300228940", // RiemannR(10^19)
            "2220819602556027015", // RiemannR(10^20)
           "21127269485932299723", // RiemannR(10^21)
          "201467286689188773625", // RiemannR(10^22)
         "1925320391607837268776", // RiemannR(10^23)
        "18435599767347541878146", // RiemannR(10^24)
       "176846309399141934626965", // RiemannR(10^25)
      "1699246750872419991992147", // RiemannR(10^26)
     "16352460426841662910939464", // RiemannR(10^27)
    "157589269275973235652219770", // RiemannR(10^28)
   "1520698109714271830281953370", // RiemannR(10^29)
  "14692398897720432716641650390"  // RiemannR(10^30)
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
  {
    for (int64_t x = 0; x < (int64_t) RiemannR_tiny.size(); x++)
    {
      std::cout << "RiemannR(" << x << ") = " << RiemannR(x);
      check(RiemannR(x) == RiemannR_tiny[x]);
    }
  }

  {
    int64_t x = 10;
    for (size_t i = 0; i < RiemannR_table.size(); i++)
    {
      std::cout << "RiemannR(" << x << ") = " << RiemannR(x);
      check(RiemannR(x) == RiemannR_table[i]);
      x *= 10;
    }
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  {
    int128_t x = ipow<15>((int128_t) 10);
    for (size_t i = 0; i < RiemannR_f128.size(); i++)
    {
      std::ostringstream oss;
      oss << RiemannR(x);
      std::cout << "RiemannR(" << x << ") = " << oss.str();
      check(oss.str() == RiemannR_f128[i]);
      x *= 10;
    }
  }

#endif

  {
    int64_t x = 10;
    for (size_t i = 0; i < RiemannR_table.size(); i++)
    {
      int64_t y = RiemannR_table[i];
      std::cout << "RiemannR_inverse(" << y << ") = " << RiemannR_inverse(y);
      check(RiemannR_inverse(y) < x &&
            RiemannR_inverse(y + 1) >= x);
      x *= 10;
    }
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  {
    int128_t x = ipow<15>((int128_t) 10);
    for (size_t i = 0; i < RiemannR_f128.size(); i++)
    {
      int128_t y = calculator::eval<int128_t>(RiemannR_f128[i]);
      std::cout << "RiemannR_inverse(" << y << ") = " << RiemannR_inverse(y);
      check(RiemannR_inverse(y) < x &&
            RiemannR_inverse(y + 1) >= x);
      x *= 10;
    }
  }

#endif

  {
    // RiemannR(8013) = 1010.00064
    int64_t x = 8013;
    int64_t y = 1010;
    std::cout << "RiemannR(" << x << ") = " << RiemannR(x);
    check(RiemannR(x) == y);

    std::cout << "RiemannR_inverse(" << y << ") = " << RiemannR_inverse(y);
    check(RiemannR_inverse(y) < x &&
          RiemannR_inverse(y + 1) >= x);
  }

  {
    // RiemannR(9557) = 1178.99908
    int64_t x = 9557;
    int64_t y = 1178;
    std::cout << "RiemannR(" << x << ") = " << RiemannR(x);
    check(RiemannR(x) == y);

    std::cout << "RiemannR_inverse(" << y << ") = " << RiemannR_inverse(y);
    check(RiemannR_inverse(y) < x &&
          RiemannR_inverse(y + 1) >= x);
  }

  // Sanity checks for tiny values of RiemannR(x)
  int64_t x;
  for (x = 0; x < 10000; x++)
  {
    int64_t rix = RiemannR(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "RiemannR(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of RiemannR(x)
  for (; x < 100000; x += 101)
  {
    int64_t rix = RiemannR(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "RiemannR(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for tiny values of RiemannR_inverse(x)
  for (x = 2; x < 1000; x++)
  {
    int64_t res = RiemannR_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "RiemannR_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of RiemannR_inverse(x)
  for (; x < 100000; x += 101)
  {
    int64_t res = RiemannR_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "RiemannR_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  {
    int64_t x = std::numeric_limits<int64_t>::max() / 10;
    int64_t res = RiemannR_inverse(x);
    if (res != std::numeric_limits<int64_t>::max())
    {
      std::cout << "RiemannR_inverse(" << x << ") != INT64_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }

#if defined(HAVE_INT128_T)
  {
    int128_t x = std::numeric_limits<int128_t>::max() / 10;
    int128_t res = RiemannR_inverse(x);
    if (res != std::numeric_limits<int128_t>::max())
    {
      std::cout << "RiemannR_inverse(" << x << ") != INT128_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
