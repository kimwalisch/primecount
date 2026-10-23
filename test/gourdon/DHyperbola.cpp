///
/// @file  DHyperbola.cpp
/// @brief Check the hyperbola terms of D and their query order in each
///        sieve segment against dense partial-sieve and Moebius
///        prefixes.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <DHyperbola.hpp>
#include <generate_primes.hpp>
#include <int128_t.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>

using namespace primecount;

struct Input
{
  int64_t xp;
  int64_t y;
  int64_t prime;
  bool easy_leaves;
};

// If x/p^3 < y, the leaves x/(p*y) < m <= y also include easy
// leaves, the last such inputs have sqrt(x/p) < p^2. If x/p^3 = y,
// they are all hard leaves and x/(p*y) = p^2.
const Array<Input, 26> inputs = {{
  { 280370, 530, 23, false },
  { 280898, 530, 23, false },
  { 280899, 530, 23, false },
  { 559681, 1058, 23, true },
  { 559682, 1058, 23, false },
  { 1414562, 1682, 29, false },
  { 560739, 1058, 23, false },
  { 560740, 1058, 23, false },
  { 1416243, 1682, 29, false },
  { 292140, 541, 23, false },
  { 807302, 899, 23, false },
  { 864017, 1200, 23, false },
  { 396017, 720, 23, false },
  { 444888, 767, 23, false },
  { 444889, 767, 23, false },
  { 444890, 767, 23, false },
  { 707280, 941, 23, false },
  { 707281, 941, 23, false },
  { 707282, 941, 23, false },
  { 921599, 1200, 23, false },
  { 921600, 1200, 23, false },
  { 279840, 530, 23, true },
  { 250000, 600, 23, true },
  { 200000, 900, 23, true },
  { 150000, 1800, 23, true },
  { 700000, 1500, 29, true }
}};

void check(bool ok, const char* message)
{
  if (!ok)
  {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

struct Reference
{
  Vector<uint32_t> lpf;
  Vector<int8_t> mu;
  Vector<int64_t> phi;
  Vector<int64_t> mertens;

  Reference(int64_t limit, int64_t prime)
  {
    auto primes = generate_primes<uint32_t>(limit);
    lpf = generate_lpf(limit, primes);
    mu = generate_moebius(limit, primes);
    phi.resize(limit + 1);
    mertens.resize(limit + 1);
    phi[0] = 0;
    mertens[0] = 0;
    for (int64_t n = 1; n <= limit; n++)
    {
      phi[n] = phi[n - 1] + (n == 1 || lpf[n] >= prime);
      mertens[n] = mertens[n - 1] + (lpf[n] > prime ? mu[n] : 0);
    }
  }
};

template <typename T, typename Prime>
void check_segments(const Input& input, int64_t size,
                    const Reference& reference, const Vector<int64_t>& expected,
                    int64_t easy, const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t p = input.prime;
  int64_t b = pi[p];
  int64_t xpy = input.xp / input.y;
  int64_t stop = expected.size() - 1;
  int64_t mertens_y = reference.mertens[input.y];
  Vector<SemiprimeEvent> events;
  SemiprimeEventSorter sorter;

  if (input.easy_leaves)
  {
    T actual = easy_hyperbola_sum(T(input.xp), p, b, input.y, mertens_y, primes, pi);
    check(actual == T(easy), "Easy hyperbola terms differ");
  }

  for (int64_t low = 0; low < stop; low += size)
  {
    int64_t high = min(low + size, stop);
    int64_t phi = reference.phi[low];
    int64_t max_m = input.y;
    uint64_t previous = 0;
    auto count = [&](uint64_t position) -> uint64_t {
      check(position < uint64_t(high - low), "Hyperbola count escaped its segment");
      check(position >= previous, "Hyperbola counts are not monotonic");
      previous = position;
      return reference.phi[low + position] - phi;
    };

    T actual = hyperbola_sum(T(input.xp), p, b, input.y, low, high, phi, mertens_y,
                             primes, pi, count, events, sorter, max_m);
    int64_t wanted = expected[high] - expected[low];
    if (actual != T(wanted) || max_m != xpy)
    {
      std::cerr << "Hyperbola segment differs: p=" << p << ", low=" << low
                << ", high=" << high << ", integer bytes=" << sizeof(T)
                << ", prime bytes=" << sizeof(Prime) << '\n';
      std::exit(1);
    }
  }
}

void check(const Input& input)
{
  int64_t p = input.prime;
  int64_t square = p * p;
  int64_t xpy = input.xp / input.y;
  int64_t sqrt_xp = isqrt(input.xp);
  int64_t xy = input.xp * p / input.y;

  // hyperbola_sum() tests the prime itself, whereas D.cpp
  // computes the easy hyperbola terms for p > (x/y)^(1/3).
  check(sqrt_xp < input.y && input.y < square * p &&
        square <= input.y &&
        input.easy_leaves == (xpy < square) &&
        input.easy_leaves == (p > iroot<3>(xy)),
        "The test must use the hyperbola method");

  // All hyperbola terms have positions < y
  int64_t stop = input.y;
  Reference reference(stop, p);
  Vector<int64_t> expected(stop + 1);
  for (int64_t i = 0; i <= stop; i++)
    expected[i] = 0;

  // Accumulate independent dense events at their exact query positions.
  // The hyperbola terms below p^2 are computed once per prime.
  int64_t easy = 0;
  auto add = [&](int64_t position, int64_t value) {
    if (position < square)
      easy += value;
    else
      expected[position + 1] += value;
  };

  add(sqrt_xp, reference.mertens[sqrt_xp] * reference.phi[sqrt_xp]);
  add(xpy, -reference.mertens[input.y] * reference.phi[xpy]);
  for (int64_t n = xpy + 1; n <= sqrt_xp; n++)
    if (reference.lpf[n] >= p)
    {
      int64_t c = reference.lpf[n] > p ? reference.mu[n] : 0;
      int64_t t = input.xp / n;
      add(t, -(reference.mertens[t] + c * reference.phi[t]));
    }

  // The hyperbola terms replace the leaves xpy < m <= y
  int64_t direct = 0;
  for (int64_t m = xpy + 1; m <= input.y; m++)
    if (reference.lpf[m] > p)
      direct -= reference.mu[m] * reference.phi[input.xp / m];
  for (int64_t i = 1; i <= stop; i++)
    expected[i] += expected[i - 1];
  check(expected.back() + easy == direct, "Hyperbola sum differs from original Moebius leaves");

  PiTable pi(input.y, 1);
  auto primes32 = pi.get_primes<uint32_t>(input.y, 1);
  Vector<int64_t> primes64(primes32.size());
  for (std::size_t i = 0; i < primes32.size(); i++)
    primes64[i] = primes32[i];
  for (int64_t size : { 240, 1680 })
  {
    check_segments<uint64_t>(input, size, reference, expected, easy, primes32, pi);
    check_segments<uint64_t>(input, size, reference, expected, easy, primes64, pi);
    #if defined(HAVE_INT128_T)
      check_segments<uint128_t>(input, size, reference, expected, easy, primes32, pi);
      check_segments<uint128_t>(input, size, reference, expected, easy, primes64, pi);
    #endif
  }
}

/// D.cpp computes easy_hyperbola_sum() for all primes
/// (x/y)^(1/3) < p <= sqrt(y), hence these primes must use
/// the hyperbola method and have easy leaves.
void check_easy_primes()
{
  for (int64_t x = 1; x <= 30000; x++)
  {
    for (int64_t y = 1; y * y <= x; y++)
    {
      for (int64_t p = iroot<3>(x / y) + 1; p * p <= y; p++)
      {
        int64_t xp = x / p;
        check(isqrt(xp) < y && p * p > y / p && xp / y < p * p,
              "Prime without easy hyperbola terms");
      }
    }
  }
}

int main()
{
  check_easy_primes();
  for (const Input& input : inputs)
    check(input);
  std::cout << "All hyperbola segment tests passed successfully!" << std::endl;
}
