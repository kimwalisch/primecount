///
/// @file  DHyperbola.cpp
/// @brief Check the hyperbola terms of D and their query order in each
///        sieve segment, and the hyperbola terms of AC's C1 formula
///        in each SegmentedPiTable segment, against dense
///        partial-sieve and Moebius prefixes.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <AC_default.hpp>
#include <DHyperbola.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <SemiprimeEventSorter.hpp>
#include <sieve/Sieve.hpp>
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
// leaves, the inputs from 279840 on have sqrt(x/p) <= p^2. If
// x/p^3 = y, they are all hard leaves and x/(p*y) = p^2.
const Array<Input, 30> inputs = {{
  { 280370, 530, 23, false },
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
  { 700000, 1500, 29, true },
  // C1 boundaries: x/(p*y) = 384 and sqrt(x/p) = 512 are segment
  // boundaries, sqrt(x/p) = 523 and x/p^3 = 509 are prime and
  // sqrt(x/p) = p^2 for x/p = p^4.
  { 230400, 600, 23, true },
  { 262144, 600, 23, true },
  { 273529, 600, 23, true },
  { 269261, 600, 23, true },
  { 279841, 600, 23, true }
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
  Vector<int8_t> omega;
  Vector<int64_t> pi;
  Vector<int64_t> phi;
  Vector<int64_t> mertens;
  Vector<int64_t> coefficients;

  Reference(int64_t limit, int64_t prime)
  {
    auto primes = generate_primes<uint32_t>(limit);
    lpf = generate_lpf(limit, primes);
    mu = generate_moebius(limit, primes);
    omega.resize(limit + 1);
    pi.resize(limit + 1);
    phi.resize(limit + 1);
    mertens.resize(limit + 1);
    coefficients.resize(limit + 1);
    omega[0] = omega[1] = 0;
    pi[0] = 0;
    phi[0] = 0;
    mertens[0] = 0;
    coefficients[0] = 0;
    for (int64_t n = 1; n <= limit; n++)
    {
      if (n > 1)
        omega[n] = omega[n / lpf[n]] + 1;
      pi[n] = pi[n - 1] + (n > 1 && lpf[n] == n);
      phi[n] = phi[n - 1] + (n == 1 || lpf[n] >= prime);
      mertens[n] = mertens[n - 1] + (lpf[n] > prime ? mu[n] : 0);
      coefficients[n] = coefficients[n - 1];
      if (lpf[n] > prime && omega[n] <= 2)
        coefficients[n] += mu[n];
      if (lpf[n] >= prime && omega[n] == 3)
        coefficients[n]++;
    }
  }

  int64_t e2(int64_t t, int64_t p) const
  {
    if (t < p * p)
      return 2 * pi[t] - 2 * pi[p] + 1;
    return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * pi[p] + 2;
  }
};

template <typename T, typename Prime>
void check_segments(const Input& input, int64_t size,
                    const Reference& reference, const Vector<int64_t>& expected,
                    const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t p = input.prime;
  int64_t b = pi[p];
  int64_t xpy = input.xp / input.y;
  int64_t stop = expected.size() - 1;
  int64_t coefficient = reference.coefficients[input.y];
  Vector<SemiprimeEvent> events;
  SemiprimeEventSorter sorter;

  for (int64_t low = 0; low < stop; low += size)
  {
    int64_t high = min(low + size, stop);
    int64_t phi = reference.phi[low];
    int64_t max_m = input.y;
    Sieve sieve(low, size, b);
    sieve.pre_sieve(primes, b - 1, low, high);
    sieve.init_counter(low, high);

    T actual = hyperbola_sum(T(input.xp), p, b, input.y, low, high, phi, coefficient,
                             primes, pi, sieve, events, sorter, max_m);
    int64_t wanted = expected[high] - expected[low];
    if (actual != T(wanted) || max_m != xpy)
    {
      std::cerr << "Hyperbola segment differs: xp=" << input.xp << ", y="
                << input.y << ", p=" << p << ", low=" << low
                << ", high=" << high << ", actual=" << actual
                << ", expected=" << T(wanted) << ", integer bytes=" << sizeof(T)
                << ", prime bytes=" << sizeof(Prime) << '\n';
      std::exit(1);
    }
  }
}

/// AC's C1 formula computes the hyperbola terms at
/// positions < p^2 in each segment [low, high[.
template <typename T, typename XP, typename Prime>
void check_c1_segments(const Input& input, int64_t size,
                       const Vector<int64_t>& easy, int64_t mertens_y,
                       const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t p = input.prime;
  int64_t b = pi[p];
  int64_t stop = easy.size() - 1;
  T x = T(input.xp) * p;
  SegmentedPiTable segmentedPi;

  for (int64_t low = 0; low < stop; low += size)
  {
    int64_t high = min(low + size, stop);
    segmentedPi.init(low, high, stop);
    T xlow = x / max(low, 1);
    T xhigh = x / high;

    T actual = C1_hyperbola(xlow, xhigh, XP(input.xp), input.y, b, mertens_y,
                            primes, pi, segmentedPi);
    int64_t wanted = easy[high] - easy[low];
    if (actual != T(wanted))
    {
      std::cerr << "C1 hyperbola segment differs: p=" << p << ", low=" << low
                << ", high=" << high << ", integer bytes=" << sizeof(T)
                << ", xp bytes=" << sizeof(XP) << ", prime bytes="
                << sizeof(Prime) << '\n';
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

  // hyperbola_sum() tests the prime itself, whereas AC computes
  // the C1 hyperbola terms for p > (x/y)^(1/3).
  check(sqrt_xp < input.y && input.y < square * square &&
        square <= input.y &&
        input.y * input.y <= input.xp * p &&
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
  // The hyperbola terms below p^2 are computed by AC's C1 formula.
  Vector<int64_t> easy(stop + 1);
  for (int64_t i = 0; i <= stop; i++)
    easy[i] = 0;
  auto add = [&](int64_t position, int64_t value) {
    if (position < square)
      easy[position + 1] += value;
    else
      expected[position + 1] += value;
  };

  add(sqrt_xp, reference.coefficients[sqrt_xp] * reference.phi[sqrt_xp]);
  add(xpy, -reference.coefficients[input.y] * reference.phi[xpy]);
  for (int64_t n = xpy + 1; n <= sqrt_xp; n++)
    if (reference.lpf[n] >= p)
    {
      int64_t c = reference.lpf[n] > p ? reference.mu[n] : 0;
      int64_t t = input.xp / n;
      add(t, reference.e2(t, p) - (1 + c) * reference.phi[t]);
    }
  // Reflect the three-factor prefix correction to positions <= H
  for (int64_t d = sqrt_xp + 1; d <= input.y; d++)
    if (reference.lpf[d] >= p && reference.omega[d] == 3)
    {
      int64_t c = reference.lpf[d] > p ? reference.mu[d] : 0;
      add(input.xp / d, (1 - c) * reference.phi[input.xp / d]);
    }

  for (int64_t t = square; t <= input.y; t++)
    check(reference.phi[t] - reference.e2(t, p) == reference.coefficients[t],
          "The dense three-factor boundary coefficient differs");

  // The hyperbola terms replace the leaves xpy < m <= y
  int64_t direct = 0;
  for (int64_t m = xpy + 1; m <= input.y; m++)
    if (reference.lpf[m] > p)
      direct -= reference.mu[m] * reference.phi[input.xp / m];
  for (int64_t i = 1; i <= stop; i++)
  {
    expected[i] += expected[i - 1];
    easy[i] += easy[i - 1];
  }
  check(expected.back() + easy.back() == direct, "Hyperbola sum differs from original Moebius leaves");

  PiTable pi(input.y, 1);
  auto primes32 = pi.get_primes<uint32_t>(input.y, 1);
  int64_t mertens_y = reference.mertens[input.y];
  Vector<int64_t> mertens = generate_mertens_y(input.y, 8, primes32, pi);
  if (input.y < square * p && pi[p] > 8)
    check(mertens[pi[p]] == mertens_y, "generate_mertens_y() differs");
  auto coefficients = generate_hyperbola_coefficients(input.y, 0, 0, primes32, pi);
  check(coefficients[pi[p]] == reference.coefficients[input.y],
        "generate_hyperbola_coefficients() differs");

  #if defined(HAVE_INT128_T)
    Vector<int64_t> primes64(primes32.size());
    for (std::size_t i = 0; i < primes32.size(); i++)
      primes64[i] = primes32[i];
    auto coefficients64 = generate_hyperbola_coefficients(input.y, 0, 0, primes64, pi);
    check(coefficients64[pi[p]] == reference.coefficients[input.y],
          "generate_hyperbola_coefficients() with 64-bit primes differs");
  #endif

  for (int64_t size : { 240, 1680 })
  {
    check_segments<uint64_t>(input, size, reference, expected, primes32, pi);
    #if defined(HAVE_INT128_T)
      check_segments<uint128_t>(input, size, reference, expected, primes32, pi);
      check_segments<uint128_t>(input, size, reference, expected, primes64, pi);
    #endif
  }

  if (input.easy_leaves)
  {
    for (int64_t size : { 128, 1024 })
    {
      check_c1_segments<uint64_t, uint64_t>(input, size, easy, mertens_y, primes32, pi);
      #if defined(HAVE_INT128_T)
        check_c1_segments<uint128_t, uint64_t>(input, size, easy, mertens_y, primes32, pi);
        check_c1_segments<uint128_t, uint128_t>(input, size, easy, mertens_y, primes64, pi);
      #endif
    }
  }
}

/// AC computes C1_hyperbola() for all primes
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
              "C1 prime without hyperbola terms below p^2");
      }
    }
  }
}

template <typename Prime>
void check_coefficients(int64_t y, const Reference& reference,
                        const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t pi_sqrty = pi[isqrt(y)];
  Vector<int64_t> mertens2(pi_sqrty + 1);
  Vector<int64_t> triples(pi_sqrty + 1);
  for (int64_t b = 1; b <= pi_sqrty; b++)
  {
    mertens2[b] = 0;
    triples[b] = 0;
    for (int64_t n = 1; n <= y; n++)
    {
      if (reference.lpf[n] > primes[b] && reference.omega[n] <= 2)
        mertens2[b] += reference.mu[n];
      if (reference.lpf[n] >= primes[b] && reference.omega[n] == 3)
        triples[b]++;
    }
  }

  const Array<int64_t, 4> k_values = { 0, 3, 8, pi_sqrty };
  const Array<int64_t, 6> min_primes = { 0, iroot<4>(y), iroot<3>(y), isqrt(y),
                                       y + 1, 2 * y };
  for (int64_t k : k_values)
    for (int64_t min_prime : min_primes)
    {
      auto coefficients = generate_hyperbola_coefficients(y, k, min_prime, primes, pi);
      check(coefficients.size() == std::size_t(pi_sqrty + 1),
            "The hyperbola coefficient vector has an incorrect size");
      for (int64_t b = k + 1; b <= pi_sqrty; b++)
      {
        int64_t expected = mertens2[b];
        if (primes[b] > max(min_prime, iroot<4>(y)))
          expected += triples[b];
        if (coefficients[b] != expected)
        {
          std::cerr << "Hyperbola coefficient differs: y=" << y << ", k=" << k
                    << ", min_prime=" << min_prime << ", p=" << primes[b]
                    << ", prime bytes=" << sizeof(Prime) << '\n';
          std::exit(1);
        }
      }
    }
}

void check_coefficients()
{
  for (int64_t y : { 80, 81, 82, 124, 125, 126, 624, 625, 626,
                     1330, 1331, 1332, 14640, 14641, 14642 })
  {
    Reference reference(y, 2);
    PiTable pi(y, 1);
    auto primes32 = pi.get_primes<uint32_t>(y, 1);
    check_coefficients(y, reference, primes32, pi);
    #if defined(HAVE_INT128_T)
      Vector<int64_t> primes64(primes32.size());
      for (std::size_t i = 0; i < primes32.size(); i++)
        primes64[i] = primes32[i];
      check_coefficients(y, reference, primes64, pi);
    #endif
  }
}

template <typename T, typename Prime>
void check_fourth_power(int64_t y, int64_t p,
                        const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t b = pi[p];
  int64_t max_m = y;
  int64_t size = 1680;
  Sieve sieve(0, size, b);
  sieve.pre_sieve(primes, b - 1, 0, size);
  sieve.init_counter(0, size);
  Vector<SemiprimeEvent> events;
  SemiprimeEventSorter sorter;
  T xp = (T(y) * y + p - 1) / p;
  T actual = hyperbola_sum(xp, p, b, y, 0, size, 0, 0,
                           primes, pi, sieve, events, sorter, max_m);
  check(actual == 0 && max_m == y,
        "The hyperbola method must reject y >= p^4");
}

void check_three_factor_inputs()
{
  // Triples on each side of H exercise all four event weights:
  // 13*17*19 and 13^3 below H, 13*17*23 and 11^2*41 above H.
  check({ 17640000, 13900, 11, false });

  for (int64_t p : { 7, 11, 23 })
    for (int64_t delta : { -1, 0, 1 })
    {
      int64_t y = p * p * p + delta;
      int64_t xp = max((y * y + p - 1) / p, y * p * p);
      check({ xp, y, p, false });
    }

  for (int64_t h : { 1330, 1331, 1332, 1680, 4198, 4199, 4200, 4320 })
    for (int64_t delta : { -1, 0, 1 })
    {
      int64_t y = h < 2000 ? 4200 : 10000;
      check({ h * h + delta, y, 11, false });
    }

  // Put repeated and square-free triples exactly at segment edges
  for (int64_t delta : { -1, 0, 1 })
  {
    check({ 997 * (1331 + delta), 2000, 11, false });
    check({ 1331 * 240 + delta, 1400, 11, false });
    check({ 2197 * 2400 + delta, 5000, 11, false });
    check({ 4199 * 4320 + delta, 10000, 11, false });
    check({ 5083 * 1680 + delta, 9000, 11, false });
    check({ 4961 * 1680 + delta, 9000, 11, false });
  }

  for (int64_t p : { 7, 11 })
  {
    int64_t fourth_power = p * p * p * p;
    int64_t y = fourth_power - 1;
    check({ (y * y + p - 1) / p, y, p, false });
    for (int64_t delta : { 0, 1 })
    {
      y = fourth_power + delta;
      PiTable pi(y, 1);
      auto primes32 = pi.get_primes<uint32_t>(y, 1);
      check_fourth_power<uint64_t>(y, p, primes32, pi);
      #if defined(HAVE_INT128_T)
        Vector<int64_t> primes64(primes32.size());
        for (std::size_t i = 0; i < primes32.size(); i++)
          primes64[i] = primes32[i];
        check_fourth_power<uint128_t>(y, p, primes32, pi);
        check_fourth_power<uint128_t>(y, p, primes64, pi);
      #endif
    }
  }
}

#if defined(HAVE_INT128_T)
template <typename Prime>
void check_large_xp(const Vector<Prime>& primes, const PiTable& pi)
{
  int64_t y = 1000;
  int64_t p = 23;
  int64_t b = pi[p];
  int64_t max_m = y;
  Sieve sieve(0, 240, b);
  sieve.pre_sieve(primes, b - 1, 0, 240);
  sieve.init_counter(0, 240);
  Vector<SemiprimeEvent> events;
  SemiprimeEventSorter sorter;
  uint128_t xp = uint128_t(1) << 100;
  uint128_t actual = hyperbola_sum(xp, p, b, y, 0, 240, 0, 0,
                                   primes, pi, sieve, events, sorter, max_m);
  check(actual == 0 && max_m == y,
        "The hyperbola method must reject a large 128-bit quotient");
}
#endif

int main()
{
  check_easy_primes();
  check_coefficients();
  for (const Input& input : inputs)
    check(input);
  check_three_factor_inputs();
  #if defined(HAVE_INT128_T)
    PiTable pi(1000, 1);
    auto primes32 = pi.get_primes<uint32_t>(1000, 1);
    auto primes64 = pi.get_primes<int64_t>(1000, 1);
    check_large_xp(primes32, pi);
    check_large_xp(primes64, pi);
  #endif
  std::cout << "All hyperbola segment tests passed successfully!" << std::endl;
}
