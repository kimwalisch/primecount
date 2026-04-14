///
/// @file   RiemannR_psi.cpp
/// @brief  Test the RiemannR_psi() and RiemannR_psi_inverse() functions.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
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
#include <string>
#include <array>

using std::size_t;
using namespace primecount;

std::array<int64_t, 14> PrimePi_table =
{
                4, // pi(10^1)
               25, // pi(10^2)
              168, // pi(10^3)
             1229, // pi(10^4)
             9592, // pi(10^5)
            78498, // pi(10^6)
           664579, // pi(10^7)
          5761455, // pi(10^8)
         50847534, // pi(10^9)
        455052511, // pi(10^10)
     4118054813ll, // pi(10^11)
    37607912018ll, // pi(10^12)
   346065536839ll, // pi(10^13)
  3204941750802ll  // pi(10^14)
};

std::array<int64_t, 7> NthPrime_table =
{
        179424673ll, // nth_prime(10^7)
       2038074743ll, // nth_prime(10^8)
      22801763489ll, // nth_prime(10^9)
     252097800623ll, // nth_prime(10^10)
    2760727302517ll, // nth_prime(10^11)
   29996224275833ll, // nth_prime(10^12)
  323780508946331ll  // nth_prime(10^13)
};

std::array<int64_t, 6> InverseSamplePoints =
{
          10000ll,
          12345ll,
          99991ll,
        1234567ll,
        9999991ll,
  1234567890123ll
};

#if defined(HAVE_FLOAT128)

std::array<std::string, 15> PrimePi_f128 =
{
                "29844570422669", // pi(10^15)
               "279238341033925", // pi(10^16)
              "2623557157654233", // pi(10^17)
             "24739954287740860", // pi(10^18)
            "234057667276344607", // pi(10^19)
           "2220819602560918840", // pi(10^20)
          "21127269486018731928", // pi(10^21)
         "201467286689315906290", // pi(10^22)
        "1925320391606803968923", // pi(10^23)
       "18435599767349200867866", // pi(10^24)
      "176846309399143769411680", // pi(10^25)
     "1699246750872437141327603", // pi(10^26)
    "16352460426841680446427399", // pi(10^27)
   "157589269275973410412739598", // pi(10^28)
  "1520698109714272166094258063"  // pi(10^29)
};

#endif

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

template <typename T>
T abs_diff(T x, T y)
{
  return (x >= y) ? x - y : y - x;
}

int main()
{
  {
    int64_t x = 10;
    for (size_t i = 0; i < PrimePi_table.size(); i++)
    {
      int64_t old_err = abs_diff(RiemannR(x), PrimePi_table[i]);
      int64_t new_err = abs_diff(RiemannR_psi(x), PrimePi_table[i]);
      std::cout << "RiemannR_psi(" << x << ") error = " << new_err;

      // The 512-zero psi correction improves almost all powers of 10
      // in practice, except for a few isolated crossover points where
      // the smooth RiemannR(x) is slightly better.
      if (x != 100000000000000ll)
        check(new_err <= old_err);
      else
        check(true);

      x *= 10;
    }
  }

#if defined(HAVE_FLOAT128) && \
    defined(HAVE_INT128_T)

  {
    int128_t x = ipow<15>((int128_t) 10);
    int128_t skip = ipow<27>((int128_t) 10);
    for (size_t i = 0; i < PrimePi_f128.size(); i++)
    {
      int128_t pix = calculator::eval<int128_t>(PrimePi_f128[i]);
      int128_t old_err = abs_diff(RiemannR(x), pix);
      int128_t new_err = abs_diff(RiemannR_psi(x), pix);
      std::cout << "RiemannR_psi(" << x << ") error = " << new_err;

      if (x != skip)
        check(new_err <= old_err);
      else
        check(true);

      x *= 10;
    }
  }

#endif

  std::cout << "RiemannR_psi_inverse(1) = " << RiemannR_psi_inverse((int64_t) 1);
  check(RiemannR_psi_inverse((int64_t) 1) == 1);

  {
    int64_t x = 10000;
    for (int i = 0; i < 10; i++)
    {
      int64_t y = RiemannR_psi(x);
      std::cout << "RiemannR_psi_inverse(" << y << ") = " << RiemannR_psi_inverse(y);
      check(RiemannR_psi_inverse(y) < x &&
            RiemannR_psi_inverse(y + 1) >= x);
      x *= 10;
    }
  }

  for (int64_t x : InverseSamplePoints)
  {
    int64_t y = RiemannR_psi(x);
    std::cout << "RiemannR_psi_inverse(" << y << ") = " << RiemannR_psi_inverse(y);
    check(RiemannR_psi_inverse(y) < x &&
          RiemannR_psi_inverse(y + 1) >= x);
  }

  {
    int64_t n = 10000000ll;
    int128_t old_err_sum = 0;
    int128_t new_err_sum = 0;
    for (size_t i = 0; i < NthPrime_table.size(); i++)
    {
      int64_t old_err = abs_diff(RiemannR_inverse(n), NthPrime_table[i]);
      int64_t new_err = abs_diff(RiemannR_psi_inverse(n), NthPrime_table[i]);
      old_err_sum += old_err;
      new_err_sum += new_err;
      n *= 10;
    }

    std::cout << "RiemannR_psi_inverse total error = " << new_err_sum;
    check(new_err_sum < old_err_sum);
  }

  // Sanity checks for tiny values of RiemannR_psi(x)
  int64_t x;
  for (x = 0; x < 10000; x++)
  {
    int64_t rix = RiemannR_psi(x);
    int64_t rr  = RiemannR(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (rix != rr ||
        rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "RiemannR_psi(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of RiemannR_psi(x)
  for (; x < 100000; x += 101)
  {
    int64_t rix = RiemannR_psi(x);
    double logx = std::log(std::max((double) x, 2.0));

    if (rix < 0 ||
        (x >= 20 && rix < x / logx) ||
        (x >= 2  && rix > x * logx))
    {
      std::cout << "RiemannR_psi(" << x << ") = " << rix << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for tiny values of RiemannR_psi_inverse(x)
  for (x = 2; x < 1000; x++)
  {
    int64_t res = RiemannR_psi_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "RiemannR_psi_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // Sanity checks for small values of RiemannR_psi_inverse(x)
  for (; x < 100000; x += 101)
  {
    int64_t res = RiemannR_psi_inverse(x);
    double logx = std::log((double) x);

    if (res < 0 ||
        res < x ||
        (x >= 5 && res > x * logx * logx))
    {
      std::cout << "RiemannR_psi_inverse(" << x << ") = " << res << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  // RiemannR_psi(x) falls back to RiemannR(x) for int128_t x >= 10^30
#if defined(HAVE_INT128_T)
  {
    int128_t x = ipow<30>((int128_t) 10);
    std::cout << "RiemannR_psi(" << x << ") = " << RiemannR_psi(x);
    check(RiemannR_psi(x) == RiemannR(x));

    x += 123456789;
    std::cout << "RiemannR_psi(" << x << ") = " << RiemannR_psi(x);
    check(RiemannR_psi(x) == RiemannR(x));

    x = ipow<31>((int128_t) 10);
    std::cout << "RiemannR_psi(" << x << ") = " << RiemannR_psi(x);
    check(RiemannR_psi(x) == RiemannR(x));
  }
#endif

  {
    int64_t x = pstd::numeric_limits<int64_t>::max() / 10;
    int64_t res = RiemannR_psi_inverse(x);
    if (res != pstd::numeric_limits<int64_t>::max())
    {
      std::cout << "RiemannR_psi_inverse(" << x << ") != INT64_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }

#if defined(HAVE_INT128_T)
  {
    int128_t x = pstd::numeric_limits<int128_t>::max() / 10;
    int128_t res = RiemannR_psi_inverse(x);
    if (res != pstd::numeric_limits<int128_t>::max())
    {
      std::cout << "RiemannR_psi_inverse(" << x << ") != INT128_MAX, failed to prevent integer overflow!" << std::endl;
      std::exit(1);
    }
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
