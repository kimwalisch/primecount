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

using std::max;
using std::size_t;
using namespace primecount;

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
       3204941731601ll, // RiemannR(10^14)
      29844570495886ll, // RiemannR(10^15)
     279238341360977ll, // RiemannR(10^16)
    2623557157055978ll, // RiemannR(10^17)
   24739954284239494ll  // RiemannR(10^18)
};

#if defined(HAVE_FLOAT128)

std::vector<std::string> RiemannR_f128 =
{
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
  int64_t x = 1;
  for (size_t i = 0; i < RiemannR_table.size(); i++)
  {
    // The accuracy of RiemannR(x) depends on
    // the width of the long double type.
    if (i >= std::numeric_limits<long double>::digits10)
      break;

    x *= 10;
    std::cout << "RiemannR(" << x << ") = " << RiemannR(x);
    check(RiemannR(x) == RiemannR_table[i]);
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  int128_t x128 = ipow((int128_t) 10, 19);
  for (size_t i = 0; i < RiemannR_f128.size(); i++)
  {
    std::ostringstream oss;
    oss << RiemannR(x128);
    std::cout << "RiemannR(" << x128 << ") = " << oss.str();
    check(oss.str() == RiemannR_f128[i]);
    x128 *= 10;
  }

#endif

  x = 1;
  for (size_t i = 0; i < RiemannR_table.size(); i++)
  {
    // The accuracy of RiemannR(x) depends on
    // the width of the long double type.
    if (i >= std::numeric_limits<long double>::digits10)
      break;

    x *= 10;
    std::cout << "RiemannR_inverse(" << RiemannR_table[i] << ") = " << RiemannR_inverse(RiemannR_table[i]);
    check(RiemannR_inverse(RiemannR_table[i]) < x &&
          RiemannR_inverse(RiemannR_table[i] + 1) >= x);
  }

  // Sanity checks for tiny values of RiemannR(x)
  for (x = 0; x < 10000; x++)
  {
    int64_t rix = RiemannR(x);
    double logx = std::log(max((double) x, 2.0));

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
    double logx = std::log(max((double) x, 2.0));

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
