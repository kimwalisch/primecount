///
/// @file  DHyperbola.cpp
/// @brief Compute the special leaves x/(p*m) with x/(p*y) < m <= y
///        of primecount's D formula using a weighted Dirichlet
///        hyperbola rearrangement, in which the sieve queries of the
///        prime leaves cancel.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
/// Copyright (C) 2026 Chopper Lin, <choppers@126.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "DHyperbola.hpp"
#include "SemiprimeEventSorter.hpp"

#include <fast_div.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <isqrt.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <sieve/Sieve.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace {

using namespace primecount;

/// The boundary coefficient is M2(y) + Q(y), where M2 is the
/// Moebius prefix restricted to at most two prime factors > p,
/// and Q counts products of three primes >= p, with multiplicity.
/// For y < p^3 this is just the usual restricted Moebius prefix.
///
template <typename Primes>
Vector<int64_t> generate_hyperbola_coefficients(int64_t y,
                                                int64_t k,
                                                int64_t min_prime,
                                                const Primes& primes,
                                                const PiTable& pi)
{
  auto coefficients = generate_mertens_y(y, k, primes, pi);
  min_prime = min(min_prime, iroot<3>(y));
  int64_t min_b = max(k, pi[max(min_prime, iroot<4>(y))]);
  int64_t triples = 0;

  for (int64_t b = pi[iroot<3>(y)]; b > min_b; b--)
  {
    int64_t yp = y / primes[b];
    int64_t max_i = pi[isqrt(yp)];

    for (int64_t i = b; i <= max_i; i++)
      triples += pi[yp / primes[i]] - i + 1;

    coefficients[b] += triples;
  }

  return coefficients;
}

/// Sum pi(xp / primes[j]) for i <= j <= last
template <typename T, typename Primes>
ALWAYS_INLINE T sum_pi(T xp,
                       int64_t i,
                       int64_t last,
                       const Primes& primes,
                       const PiTable& pi)
{
  T sum = 0;

  // Unroll loop to increase instruction level parallelism
  for (; i + 3 <= last; i += 4)
  {
    int64_t xpq0 = fast_div64(xp, primes[i]);
    int64_t xpq1 = fast_div64(xp, primes[i + 1]);
    int64_t xpq2 = fast_div64(xp, primes[i + 2]);
    int64_t xpq3 = fast_div64(xp, primes[i + 3]);
    sum += pi[xpq0] + pi[xpq1] + pi[xpq2] + pi[xpq3];
  }

  NO_UNROLL_LOOP
  for (; i <= last; i++)
    sum += pi[fast_div64(xp, primes[i])];

  return sum;
}

/// Returns Phi2(t) - M2(t) for p = primes[b] and t >= p^2,
/// restricting both prefixes to at most two prime factors.
/// This equals phi(t, b - 1) - M(t) when t < p^3.
///
ALWAYS_INLINE int64_t phi_minus_mertens(int64_t t,
                                        int64_t p,
                                        int64_t b,
                                        const PiTable& pi)
{
  return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * b + 2;
}

/// Sum phi_minus_mertens(X/q) over the primes q in (min_n, max_n],
/// where X = x/p and all X/q are >= p^2. Gourdon's inversion
/// equality turns the sums of pi(X/(p*q)) and pi(sqrt(X/q)) into sums
/// over the reflected primes r, avoiding two pi lookups per prime q.
/// Use INVERT = false for the small third primes of triples.
///
template <bool INVERT = true, typename T, typename Primes>
T sum_phi_minus_mertens(T xp,
                        int64_t p,
                        int64_t b,
                        int64_t min_n,
                        int64_t max_n,
                        const Primes& primes,
                        const PiTable& pi)
{
  int64_t pi_min_n = pi[min_n];
  int64_t pi_max_n = pi[max_n];

  if (pi_min_n == pi_max_n)
    return 0;

  if (!INVERT)
  {
    T sum = 0;
    int64_t i = pi_min_n + 1;

    // The third primes of triples are small, so reflecting their
    // interval would expand it. Evaluate their E2 terms directly.
    for (; i + 3 <= pi_max_n; i += 4)
    {
      int64_t t0 = fast_div64(xp, primes[i]);
      int64_t t1 = fast_div64(xp, primes[i + 1]);
      int64_t t2 = fast_div64(xp, primes[i + 2]);
      int64_t t3 = fast_div64(xp, primes[i + 3]);
      sum += phi_minus_mertens(t0, p, b, pi);
      sum += phi_minus_mertens(t1, p, b, pi);
      sum += phi_minus_mertens(t2, p, b, pi);
      sum += phi_minus_mertens(t3, p, b, pi);
    }

    NO_UNROLL_LOOP
    for (; i <= pi_max_n; i++)
      sum += phi_minus_mertens(fast_div64(xp, primes[i]), p, b, pi);

    return sum;
  }

  T xpp = xp / p;
  int64_t min_r = fast_div64(xpp, max_n);
  int64_t max_r = fast_div64(xpp, min_n + 1);
  int64_t min_r2 = isqrt(fast_div64(xp, max_n));
  int64_t max_r2 = isqrt(fast_div64(xp, min_n + 1));

  int64_t pi_min_r = pi[min_r];
  int64_t pi_max_r = pi[max_r];
  int64_t pi_min_r2 = pi[min_r2];
  int64_t pi_max_r2 = pi[max_r2];

  T sum = sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
  sum += sum_pi(xpp, pi_min_r + 1, pi_max_r, primes, pi);
  sum += T(pi_max_n) * pi_min_r - T(pi_min_n) * pi_max_r;
  sum += T(pi_max_n) * pi_min_r2 - T(pi_min_n) * pi_max_r2;
  sum -= T(pi_max_n - pi_min_n) * (4 * b - 2);

  for (int64_t j = pi_min_r2 + 1; j <= pi_max_r2; j++)
  {
    int64_t r = primes[j];
    sum += pi[fast_div64(xp, r * r)];
  }

  return sum;
}

/// Generate triples q*r*s in (min_n, max_n], p <= q <= r <= s.
/// Above sqrt(X) their reflected correction has weight a-c.
/// Below sqrt(X) all triples contribute E2(X/n), but only those
/// containing p or a repeated prime need a sieve query (weight 1).
///
template <bool UPPER, typename T, typename Primes>
T triple_events(T xp,
                int64_t b,
                int64_t min_n,
                int64_t max_n,
                const Primes& primes,
                const PiTable& pi,
                Vector<SemiprimeEvent>& events)
{
  T sum = 0;
  int64_t max_i = pi[iroot<3>(max_n)];
  int64_t p = primes[b];

  for (int64_t i = b; i <= max_i; i++)
  {
    int64_t q = primes[i];
    int64_t max_j = pi[isqrt(max_n / q)];

    for (int64_t j = i; j <= max_j; j++)
    {
      int64_t qr = q * primes[j];
      int64_t min_k = max(j, pi[min_n / qr] + 1);
      int64_t max_k = pi[max_n / qr];

      if (!UPPER && i != b && i != j)
      {
        int64_t min_s = max(min_n / qr, int64_t(primes[j]));
        int64_t max_s = max_n / qr;

        if (max(min_k, j + 1) <= max_k)
          sum += sum_phi_minus_mertens<false>(xp / qr, p, b, min_s, max_s, primes, pi);

        // Only repeated factors or a factor p have a+c != 0
        max_k = min(max_k, j);
      }

      for (int64_t k = max_k; k >= min_k; k--)
      {
        int64_t n = qr * primes[k];

        if (UPPER)
        {
          int64_t weight = (i != b && i != j && j != k) + 1;
          events.push_back(SemiprimeEvent{n, weight});
        }
        else
          events.push_back(SemiprimeEvent{n, 1});
      }
    }
  }

  return sum;
}

/// Compute the special leaves x/(p*y) < m <= y of the prime
/// p = primes[b] in the current segment [low, high[ using the
/// hyperbola method, if sqrt(x/p) < y < p^4. Then max_m is lowered
/// to x/(p*y) and the caller processes the remaining leaves using
/// the factor table.
///
/// Let X = floor(x/p), L = xpy = floor(X/y), H = sqrt_xp =
/// floor(sqrt(X)) and Phi(t) = phi(t, b - 1). The leaves L < m <= y
/// contribute
///   (Phi(H) - E2(H)) Phi(H) - (M2(y) + Q(y)) Phi(L)
///   + sum_{L < n <= H} a(n) E2(X/n)
///   - sum_{L < n <= H} (a(n) + c(n)) Phi(X/n)
///   + sum_{H < n <= y, Omega(n) = 3} (a(n) - c(n)) Phi(X/n),
/// where a(n) = [lpf(n) >= p], c(n) = mu(n)[lpf(n) > p] and
/// E2 = Phi2 - M2 and Q counts triples with prime factors >= p.
/// Reflecting the three-factor part of E = Phi - M gives the last
/// sum and changes the two boundary coefficients. This avoids
/// evaluating a three-factor prefix for every cancelled prime.
/// Only terms whose position (Phi argument) is in
/// [max(low, p^2), high[ are computed here, the terms at positions
/// < p^2 are computed by the C1 formula, see C1_hyperbola() in
/// AC_default.hpp.
/// Prime n and square-free triples with lpf(n) > p have a+c = 0.
/// Other composites need a query with weight 1 or 2. Each list of
/// events is processed in decreasing n, as count() requires
/// increasing positions X/n. See doc/Hard-Special-Leaves-Hyperbola.md.
///
template <typename T, typename Primes>
ALWAYS_INLINE T hyperbola_sum(T xp,
                              int64_t p,
                              int64_t b,
                              int64_t y,
                              int64_t low,
                              int64_t high,
                              int64_t phi,
                              int64_t coefficient,
                              const Primes& primes,
                              const PiTable& pi,
                              Sieve& sieve,
                              Vector<SemiprimeEvent>& events,
                              SemiprimeEventSorter& sorter,
                              int64_t& max_m)
{
  int64_t sqrt_xp = isqrt(xp);

  if (sqrt_xp >= y)
    return 0;

  int64_t square = p * p;
  bool triples = square <= y / p;

  if (triples && square <= y / square)
    return 0;

  int64_t xpy = fast_div64(xp, y);

  // All positions are >= L. Using it here also keeps the
  // quotient X/query_low below 2*y, even in the 128-bit case.
  int64_t query_low = max3(low, square, xpy);
  max_m = min(max_m, xpy);
  T sum = 0;

  if (high <= xpy)
    return 0;

  if (query_low <= xpy && xpy < high)
  {
    int64_t phi_xpy = phi + sieve.count(xpy - low);
    sum -= T(phi_xpy) * coefficient;
  }

  // The reflected three-factor terms have positions <= H,
  // so their queries must precede the boundary query at H.
  if (triples && query_low <= sqrt_xp)
  {
    int64_t min_n = min(fast_div(xp, high), y);
    int64_t max_n = min(fast_div64(xp, query_low), y);
    min_n = max(min_n, sqrt_xp);

    if (min_n < max_n)
    {
      events.clear();
      triple_events<true>(xp, b, min_n, max_n, primes, pi, events);
      sorter.sort(events, min_n, max_n);

      for (const SemiprimeEvent& event : events)
      {
        int64_t xpn = fast_div64(xp, event.n);
        sum += event.weight * (phi + sieve.count(xpn - low));
      }
    }
  }

  if (query_low <= sqrt_xp && sqrt_xp < high)
  {
    int64_t phi_sqrt_xp = phi + sieve.count(sqrt_xp - low);
    sum += T(phi_sqrt_xp) * (phi_sqrt_xp - phi_minus_mertens(sqrt_xp, p, b, pi));
  }

  // n in (L, H] with position X/n in [query_low, high[
  int64_t min_n = min(fast_div(xp, high), sqrt_xp);
  int64_t max_n = min(fast_div64(xp, query_low), sqrt_xp);
  min_n = max(min_n, xpy);

  if (min_n < max_n)
  {
    sum += sum_phi_minus_mertens(xp, p, b, min_n, max_n, primes, pi);
    int64_t max_i = pi[isqrt(max_n)];
    events.clear();

    for (int64_t i = b; i <= max_i; i++)
    {
      int64_t q = primes[i];
      int64_t min_j = max(i, pi[min_n / q] + 1);
      int64_t max_j = pi[max_n / q];

      for (int64_t j = max_j; j >= min_j; j--)
      {
        int64_t n = q * primes[j];
        int64_t weight = (i != b && i != j) + 1;
        events.push_back(SemiprimeEvent{n, weight});
      }
    }

    if (square <= max_n / p)
      sum += triple_events<false>(xp, b, min_n, max_n, primes, pi, events);

    if (events.empty())
      return sum;

    sorter.sort(events, min_n, max_n);

    // The positions t = X/n increase, hence the
    // indexes pi(t/p) and pi(sqrt(t)) never decrease.
    int64_t t = fast_div64(xp, events[0].n);
    int64_t pi_tp = pi[t / p];
    int64_t pi_sqrt = pi[isqrt(t)];
    int64_t next_tp = primes[pi_tp + 1] * p;
    int64_t r = primes[pi_sqrt + 1];
    int64_t next_sqrt = r * r;

    for (const SemiprimeEvent& event : events)
    {
      int64_t xpn = fast_div64(xp, event.n);
      int64_t phi_xpn = phi + sieve.count(xpn - low);
      sum -= event.weight * phi_xpn;

      // Compute phi_minus_mertens(xpn) by incrementally updating
      // pi(xpn/p) and pi(sqrt(xpn)) instead of computing them from
      // scratch. Summed over all segments, the loops iterate at
      // most pi(y/p) + pi(sqrt(y)) times for the prime p, and over
      // all primes this is amortized O(1) per event. Hence the
      // runtime complexity is unchanged, but most of the time they
      // don't iterate at all, which saves operations.

      NO_UNROLL_LOOP
      NO_VECTORIZE_LOOP
      while (xpn >= next_tp)
      {
        pi_tp++;
        next_tp = primes[pi_tp + 1] * p;
      }

      NO_UNROLL_LOOP
      NO_VECTORIZE_LOOP
      while (xpn >= next_sqrt)
      {
        pi_sqrt++;
        r = primes[pi_sqrt + 1];
        next_sqrt = r * r;
      }

      // phi_minus_mertens(xpn), its -4b + 2 term is added below
      sum += 2 * pi[xpn] + pi_tp + pi_sqrt;
    }

    sum -= events.size() * T(4 * b - 2);
  }

  return sum;
}

} // namespace

namespace primecount {

Vector<int64_t> generate_hyperbola_coefficients(int64_t y,
                                                int64_t k,
                                                int64_t min_prime,
                                                const Vector<uint32_t>& primes,
                                                const PiTable& pi)
{
  return ::generate_hyperbola_coefficients(y, k, min_prime, primes, pi);
}

Vector<int64_t> generate_hyperbola_coefficients(int64_t y,
                                                int64_t k,
                                                int64_t min_prime,
                                                const Vector<int64_t>& primes,
                                                const PiTable& pi)
{
  return ::generate_hyperbola_coefficients(y, k, min_prime, primes, pi);
}

uint64_t hyperbola_sum(uint64_t xp,
                       int64_t p,
                       int64_t b,
                       int64_t y,
                       int64_t low,
                       int64_t high,
                       int64_t phi,
                       int64_t coefficient,
                       const Vector<uint32_t>& primes,
                       const PiTable& pi,
                       Sieve& sieve,
                       Vector<SemiprimeEvent>& events,
                       SemiprimeEventSorter& sorter,
                       int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, coefficient,
                         primes, pi, sieve, events, sorter, max_m);
}

#if defined(HAVE_INT128_T)

uint128_t hyperbola_sum(uint128_t xp,
                        int64_t p,
                        int64_t b,
                        int64_t y,
                        int64_t low,
                        int64_t high,
                        int64_t phi,
                        int64_t coefficient,
                        const Vector<uint32_t>& primes,
                        const PiTable& pi,
                        Sieve& sieve,
                        Vector<SemiprimeEvent>& events,
                        SemiprimeEventSorter& sorter,
                        int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, coefficient,
                         primes, pi, sieve, events, sorter, max_m);
}

uint128_t hyperbola_sum(uint128_t xp,
                        int64_t p,
                        int64_t b,
                        int64_t y,
                        int64_t low,
                        int64_t high,
                        int64_t phi,
                        int64_t coefficient,
                        const Vector<int64_t>& primes,
                        const PiTable& pi,
                        Sieve& sieve,
                        Vector<SemiprimeEvent>& events,
                        SemiprimeEventSorter& sorter,
                        int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, coefficient,
                         primes, pi, sieve, events, sorter, max_m);
}

#endif

} // namespace
