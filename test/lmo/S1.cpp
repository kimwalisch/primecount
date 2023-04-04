///
/// @file   S1.cpp
/// @brief  Test the computation of the ordinary leaves
///         S1(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <imath.hpp>
#include <S.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, 10000000);
  int threads = 1;

  for (int i = 0; i < 1000; i++)
  {
    int64_t x = dist(gen);
    int64_t y = iroot<3>(x);
    int64_t c = PhiTiny::get_c(y);
    int64_t s1 = 0;

    auto primes = generate_n_primes<int32_t>(c);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // ordinary leaves
    for (int64_t n = 1; n <= y; n++)
      if (lpf[n] > primes[c])
        s1 += mu[n] * phi_tiny(x / n, c);

    std::cout << "S1(" << x << ", " << y << ") = " << s1;
    check(s1 == S1(x, y, c, threads));
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
