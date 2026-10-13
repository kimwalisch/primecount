///
/// @file  DCentral.cpp
/// @brief Check central D contributions and query order in each sieve
///        segment against dense partial-sieve and Moebius prefixes.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <DCentral.hpp>
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
  int64_t z;
  int64_t prime;
  bool natural;
};

const Array<Input, 22> inputs = {{
  { 280370, 530, 530, 23, true },
  { 280898, 530, 530, 23, true },
  { 280899, 530, 530, 23, false },
  { 559681, 1058, 1058, 23, true },
  { 559682, 1058, 1058, 23, true },
  { 280370, 600, 700, 23, true },
  { 1414562, 1682, 1682, 29, true },
  { 560739, 1058, 1058, 23, false },
  { 560740, 1058, 1058, 23, false },
  { 1416243, 1682, 1682, 29, false },
  { 292140, 541, 541, 23, false },
  { 807302, 899, 899, 23, false },
  { 864017, 1200, 1200, 23, false },
  { 396017, 720, 720, 23, false },
  { 444888, 767, 767, 23, false },
  { 444889, 767, 767, 23, false },
  { 444890, 767, 767, 23, false },
  { 707280, 941, 941, 23, false },
  { 707281, 941, 941, 23, false },
  { 707282, 941, 941, 23, false },
  { 921599, 1200, 1200, 23, false },
  { 921600, 1200, 1200, 23, false }
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
                    const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t p = input.prime;
  int64_t b = pi[p];
  int64_t upper = input.natural ? input.xp / (p * p) : input.y;
  int64_t lower = input.xp / upper;
  int64_t middle = isqrt(input.xp);
  int64_t stop = expected.size() - 1;
  Vector<CentralEvent> events;
  SemiprimeEventBitmap sorter;

  for (int64_t low = 0; low < stop; low += size)
  {
    int64_t high = min(low + size, stop);
    int64_t phi = reference.phi[low];
    uint64_t previous = 0;
    auto count = [&](uint64_t position) -> uint64_t {
      check(position < uint64_t(high - low), "Central count escaped its segment");
      check(position >= previous, "Central counts are not monotonic");
      previous = position;
      return reference.phi[low + position] - phi;
    };

    T actual;
    if (input.natural)
      actual = central_square_segment(T(input.xp), p, b, upper, middle,
                                      input.z, low, high, phi, primes, pi,
                                      count, events, sorter);
    else
      actual = central_y_segment(T(input.xp), p, b, lower, middle,
                                 reference.mertens[upper], low, high, phi,
                                 primes, pi, count, events, sorter);
    int64_t wanted = expected[high] - expected[low];
    if (actual != T(wanted))
    {
      std::cerr << "Central segment differs: p=" << p << ", low=" << low
                << ", high=" << high << ", integer bytes=" << sizeof(T)
                << ", prime bytes=" << sizeof(Prime) << '\n';
      std::exit(1);
    }
  }
}

void check(const Input& input)
{
  int64_t p = input.prime;
  int64_t upper = input.natural ? input.xp / (p * p) : input.y;
  int64_t lower = input.xp / upper;
  int64_t middle = isqrt(input.xp);
  check(p * p <= lower && lower <= middle && middle < upper &&
        upper < p * p * p && upper <= input.y,
        "The test must enter a valid central region");
  check(input.natural ? lower == p * p : input.z == input.y,
        "Central test region does not match its parameters");
  int64_t min_m = input.natural ? max(p, input.z / p) : lower;
  int64_t stop = input.xp / (min_m + 1) + 1;
  Reference reference(max(stop, input.y), p);
  Vector<int64_t> expected(stop + 1);
  for (int64_t i = 0; i <= stop; i++)
    expected[i] = 0;

  // Accumulate independent dense events at their exact query positions
  expected[middle + 1] += reference.mertens[middle] * reference.phi[middle];
  int64_t boundary = input.natural ? upper : lower;
  expected[boundary + 1] -= reference.mertens[upper] * reference.phi[lower];
  for (int64_t n = lower + 1; n <= middle; n++)
    if (reference.lpf[n] >= p)
    {
      int64_t c = reference.lpf[n] > p ? reference.mu[n] : 0;
      int64_t t = input.xp / n;
      expected[t + 1] -= reference.mertens[t] + c * reference.phi[t];
    }
  if (input.natural)
    for (int64_t m = min_m + 1; m <= lower; m++)
      if (reference.lpf[m] > p)
        expected[input.xp / m + 1] -= reference.mu[m] * reference.phi[input.xp / m];

  int64_t direct = 0;
  for (int64_t m = min_m + 1; m <= upper; m++)
    if (reference.lpf[m] > p)
      direct -= reference.mu[m] * reference.phi[input.xp / m];
  for (int64_t i = 1; i <= stop; i++)
    expected[i] += expected[i - 1];
  check(expected.back() == direct, "Central sum differs from original Moebius leaves");

  PiTable pi(input.y, 1);
  auto primes32 = pi.get_primes<uint32_t>(input.y, 1);
  Vector<int64_t> primes64(primes32.size());
  for (std::size_t i = 0; i < primes32.size(); i++)
    primes64[i] = primes32[i];
  for (int64_t size : { 240, 1680 })
  {
    check_segments<uint64_t>(input, size, reference, expected, primes32, pi);
    check_segments<uint64_t>(input, size, reference, expected, primes64, pi);
    #if defined(HAVE_INT128_T)
      check_segments<uint128_t>(input, size, reference, expected, primes32, pi);
      check_segments<uint128_t>(input, size, reference, expected, primes64, pi);
    #endif
  }
}

int main()
{
  for (const Input& input : inputs)
    check(input);
  std::cout << "All central-region segment tests passed successfully!" << std::endl;
}
