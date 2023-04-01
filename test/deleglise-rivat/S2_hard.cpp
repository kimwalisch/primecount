///
/// @file   S2_hard.cpp
/// @brief  Test the computation of the hard special leaves
///         S2_hard(x, y) used in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime counting algorithms.
///
///         Note: when we set y = x^(1/3) then there are no
///         trivial and no easy special leaves which allows
///         us to test only the hard special leaves.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
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
  std::uniform_int_distribution<int> dist(1, 10000000);
  int threads = 1;

  // test small x
  for (int i = 1; i < 30000; i++)
  {
    int64_t x = i;
    int64_t y = iroot<3>(x);
    int64_t pi_y = pi_noprint(y, threads);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    std::cout << "S2_hard(" << x << ", " << y << ") = " << s2;
    check(s2 == S2_hard(x, y, z, c, Li(x), threads));
  }

  // test random x
  for (int i = 0; i < 500; i++)
  {
    int64_t x = dist(gen);
    int64_t y = iroot<3>(x);
    int64_t pi_y = pi_noprint(y, threads);
    int64_t z = x / y;
    int64_t c = PhiTiny::get_c(y);
    int64_t s2 = 0;

    auto primes = generate_primes<int32_t>(y);
    auto lpf = generate_lpf(y);
    auto mu = generate_moebius(y);

    // special leaves
    for (int64_t b = c + 1; b < pi_y; b++)
      for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
        if (lpf[m] > primes[b])
          s2 -= mu[m] * phi(x / (primes[b] * m), b - 1);

    std::cout << "S2_hard(" << x << ", " << y << ") = " << s2;
    check(s2 == S2_hard(x, y, z, c, Li(x), threads));
  }

  threads = get_num_threads();

  {
    // Test S2_hard(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 178815;
    int64_t z = 55923720;
    int64_t c = 8;
    int64_t res1 = S2_hard(x, y, z, c, Li(x), threads);
    int64_t res2 = 371595005834ll;

    std::cout << "S2_hard(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test S2_hard(1e14) and compare with known correct value
    int64_t x = 100000000000000ll;
    int64_t y = 494134;
    int64_t z = 202374254;
    int64_t c = 8;
    int64_t res1 = S2_hard(x, y, z, c, Li(x), threads);
    int64_t res2 = 3474606376629ll;

    std::cout << "S2_hard(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test S2_hard(1e15) and compare with known correct value
    int128_t x = 1000000000000000ll;
    int64_t y = 1378500;
    int64_t z = 725426187;
    int64_t c = 8;
    int128_t res1 = S2_hard(x, y, z, c, Li(x), threads);
    int128_t res2 = 32272468597257ll;

    std::cout << "S2_hard(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
