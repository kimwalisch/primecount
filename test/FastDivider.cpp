///
/// @file  FastDivider.cpp
/// @brief Test FastDivider, i.e. the fast division of numerators
///        < 2^63 by multiplication with precomputed reciprocals.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <random>

#if defined(ENABLE_FASTDIVIDER)
  #include <FastDivider.hpp>
#endif

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

#if defined(ENABLE_FASTDIVIDER)

bool check_div(uint64_t n, uint64_t d, const FastDivider& fd)
{
  if (n > FastDivider::max() || n / fd == n / d)
    return true;

  std::cout << "\n" << n << " / FastDivider(" << d << ") = " << n / fd;
  std::cout << ", correct: " << n / d << "\n";
  return false;
}

/// Test the numerators for which the multiplication
/// by the reciprocal is most likely to fail.
///
bool test_divisor(uint64_t d, std::mt19937_64& gen)
{
  FastDivider fd(d);
  uint64_t max_n = FastDivider::max();
  uint64_t max_q = max_n / d;

  // Largest multiples of d and their neighbours
  for (uint64_t i = 0; i < 3 && i <= max_q; i++)
  {
    uint64_t n = (max_q - i) * d;
    if (!check_div(n, d, fd) ||
        !check_div(n + 1, d, fd))
      return false;
    if (n > 0 &&
        !check_div(n - 1, d, fd))
      return false;
    if (d - 1 <= max_n - n &&
        !check_div(n + d - 1, d, fd))
      return false;
  }

  uint64_t special[] = { 0, 1, d - 1, d, d + 1, max_n - 1, max_n,
                         (1ull << 32) - 1, 1ull << 32, (1ull << 32) + 1 };

  for (uint64_t n : special)
    if (!check_div(n, d, fd))
      return false;

  for (int i = 0; i < 8; i++)
  {
    uint64_t n1 = gen() >> 1;
    uint64_t n2 = gen() >> (1 + gen() % 63);
    uint64_t n3 = (gen() % (max_q + 1)) * d;
    if (!check_div(n1, d, fd) ||
        !check_div(n2, d, fd) ||
        !check_div(n3, d, fd))
      return false;
    if (d - 1 <= max_n - n3 &&
        !check_div(n3 + d - 1, d, fd))
      return false;
  }

  return true;
}

#endif

int main()
{
#if defined(ENABLE_FASTDIVIDER)

  std::random_device rd;
  std::mt19937_64 gen(rd());

  std::cout << "sizeof(FastDivider) = " << sizeof(FastDivider);
  check(sizeof(FastDivider) == 9);

  std::cout << "FastDivider::max() = " << FastDivider::max();
  check(FastDivider::max() == (1ull << 63) - 1);

  bool OK = true;

#if defined(HAVE_INT128_T)
  // Compare with 128-bit integer arithmetic
  for (int i = 0; i < 100000 && OK; i++)
  {
    uint64_t x = gen() >> (gen() % 64);
    uint64_t y = gen() >> (gen() % 64);
    OK = mulhi64(x, y) == uint64_t((uint128_t(x) * y) >> 64);

    if (OK && y > 0)
    {
      uint64_t numhi = x % y;
      uint64_t numlo = gen();
      uint128_t numer = (uint128_t(numhi) << 64) | numlo;
      OK = div128_64(numhi, numlo, y) == uint64_t(numer / y);
    }
  }

  std::cout << "mulhi64(x, y) and div128_64(numhi, numlo, den)";
  check(OK);
#endif

  for (uint64_t d = 2; d <= 20000 && OK; d++)
    OK = test_divisor(d, gen);

  std::cout << "n / FastDivider(d) for 2 <= d <= 20000";
  check(OK);

  // Divisors near powers of 2 and the largest divisors
  for (int k = 1; k < 64 && OK; k++)
    for (uint64_t i = 0; i <= 32 && OK; i++)
      if ((1ull << k) + i >= 18)
        OK = test_divisor((1ull << k) + i - 16, gen);

  for (uint64_t i = 0; i <= 16 && OK; i++)
    OK = test_divisor(~0ull - i, gen);

  std::cout << "n / FastDivider(d) for d near 2^k";
  check(OK);

  // Random divisors of all bit sizes
  for (int bits = 2; bits <= 64 && OK; bits++)
  {
    for (int i = 0; i < 1000 && OK; i++)
    {
      uint64_t d = gen() >> (64 - bits);
      d |= 1ull << (bits - 1);
      OK = test_divisor(d, gen);
    }
  }

  std::cout << "n / FastDivider(d) for random d";
  check(OK);

#else
  std::cout << "FastDivider is disabled, skipping tests" << std::endl;
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
