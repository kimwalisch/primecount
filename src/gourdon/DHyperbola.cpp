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

/// Returns E_p(t) = Phi_p(t) - M_p(t) for p = primes[b] and
/// p^2 <= t < p^3, where Phi_p(t) = phi(t, b - 1) and
/// M_p(t) = sum_{n <= t, lpf(n) > p} mu(n).
///
ALWAYS_INLINE int64_t phi_minus_mertens(int64_t t,
                                        int64_t p,
                                        int64_t b,
                                        const PiTable& pi)
{
  return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * b + 2;
}

/// Sum E_p(X/n) over the primes n in (min_n, max_n], where X = x/p
/// and all X/n are in [p^2, p^3[. Gourdon's inversion equality
/// turns the sums of pi(X/(p*n)) and pi(sqrt(X/n)) into sums over
/// the reflected primes r with p*n*r <= X and q with n*q^2 <= X,
/// avoiding two pi lookups per prime n.
///
template <typename T, typename Primes>
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

  T xpp = xp / p;
  int64_t min_r = fast_div64(xpp, max_n);
  int64_t max_r = fast_div64(xpp, min_n + 1);
  int64_t min_q = isqrt(fast_div64(xp, max_n));
  int64_t max_q = isqrt(fast_div64(xp, min_n + 1));

  int64_t pi_min_r = pi[min_r];
  int64_t pi_max_r = pi[max_r];
  int64_t pi_min_q = pi[min_q];
  int64_t pi_max_q = pi[max_q];

  T sum = sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
  sum += sum_pi(xpp, pi_min_r + 1, pi_max_r, primes, pi);
  sum += T(pi_max_n) * pi_min_r - T(pi_min_n) * pi_max_r;
  sum += T(pi_max_n) * pi_min_q - T(pi_min_n) * pi_max_q;
  sum -= T(pi_max_n - pi_min_n) * (4 * b - 2);

  for (int64_t j = pi_min_q + 1; j <= pi_max_q; j++)
  {
    int64_t q = primes[j];
    sum += pi[fast_div64(xp, q * q)];
  }

  return sum;
}

/// Compute the hyperbola terms of the prime p = primes[b] with
/// positions in [max(low, p^2), high[, if sqrt(x/p) < y < p^3.
/// Then max_m is lowered to x/(p*y) and the caller processes the
/// remaining leaves using the factor table. C1_hyperbola() in
/// AC_default.hpp computes the terms at positions < p^2.
///
/// Let X = floor(x/p), L = xpy = floor(X/y), H = sqrt_xp =
/// floor(sqrt(X)), with Phi_p, M_p and E_p as in phi_minus_mertens().
/// The leaves L < m <= y contribute the hyperbola terms
///   M_p(H) Phi_p(H) - M_p(y) Phi_p(L)
///   + sum_{L < n <= H} a_p(n) E_p(X/n)
///   - sum_{L < n <= H} (a_p(n) + c_p(n)) Phi_p(X/n),
/// where a_p(n) = [lpf(n) >= p] and c_p(n) = mu(n)[lpf(n) > p].
/// Prime n have a_p(n) + c_p(n) = 0, so their sieve queries cancel.
/// The other n (events) are q*r with p <= q <= r and weight
/// a_p(n) + c_p(n) = 1 if q = p or q = r, else 2. They are processed
/// in decreasing n, i.e. increasing positions X/n, as required by
/// sieve.count() and the incremental E_p computation.
///
template <typename T, typename Primes>
ALWAYS_INLINE T D_hyperbola(T xp,
                            int64_t p,
                            int64_t b,
                            int64_t y,
                            int64_t low,
                            int64_t high,
                            int64_t phi,
                            int64_t mertens_y,
                            const Primes& primes,
                            const PiTable& pi,
                            Sieve& sieve,
                            Vector<SemiprimeEvent>& events,
                            SemiprimeEventSorter& sorter,
                            int64_t& max_m)
{
  int64_t square = p * p;
  int64_t sqrt_xp = isqrt(xp);

  if (sqrt_xp >= y || square <= y / p)
    return 0;

  int64_t xpy = fast_div64(xp, y);
  int64_t query_low = max(low, square);
  max_m = min(max_m, xpy);
  T sum = 0;

  if (query_low <= xpy && xpy < high)
  {
    int64_t phi_xpy = phi + sieve.count(xpy - low);
    sum -= T(phi_xpy) * mertens_y;
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
      t = fast_div64(xp, event.n);
      int64_t phi_t = phi + sieve.count(t - low);
      sum -= event.weight * phi_t;

      // Compute phi_minus_mertens(t) by incrementally updating
      // pi(t/p) and pi(sqrt(t)) instead of computing them from
      // scratch. Summed over all segments, the loops iterate at
      // most pi(y/p) + pi(sqrt(y)) times for the prime p, and
      // over all primes this is amortized O(1) per event. Hence
      // the runtime complexity is unchanged, but most of the
      // time they don't iterate at all, which saves operations.

      NO_UNROLL_LOOP
      NO_VECTORIZE_LOOP
      while (t >= next_tp)
      {
        pi_tp++;
        next_tp = primes[pi_tp + 1] * p;
      }

      NO_UNROLL_LOOP
      NO_VECTORIZE_LOOP
      while (t >= next_sqrt)
      {
        pi_sqrt++;
        r = primes[pi_sqrt + 1];
        next_sqrt = r * r;
      }

      // phi_minus_mertens(t), its -4b + 2 term is added below
      sum += 2 * pi[t] + pi_tp + pi_sqrt;
    }

    sum -= events.size() * T(4 * b - 2);
  }

  return sum;
}

} // namespace

namespace primecount {

uint64_t D_hyperbola(uint64_t xp,
                     int64_t p,
                     int64_t b,
                     int64_t y,
                     int64_t low,
                     int64_t high,
                     int64_t phi,
                     int64_t mertens_y,
                     const Vector<uint32_t>& primes,
                     const PiTable& pi,
                     Sieve& sieve,
                     Vector<SemiprimeEvent>& events,
                     SemiprimeEventSorter& sorter,
                     int64_t& max_m)
{
  return ::D_hyperbola(xp, p, b, y, low, high, phi, mertens_y,
                       primes, pi, sieve, events, sorter, max_m);
}

#if defined(HAVE_INT128_T)

uint128_t D_hyperbola(uint128_t xp,
                      int64_t p,
                      int64_t b,
                      int64_t y,
                      int64_t low,
                      int64_t high,
                      int64_t phi,
                      int64_t mertens_y,
                      const Vector<uint32_t>& primes,
                      const PiTable& pi,
                      Sieve& sieve,
                      Vector<SemiprimeEvent>& events,
                      SemiprimeEventSorter& sorter,
                      int64_t& max_m)
{
  return ::D_hyperbola(xp, p, b, y, low, high, phi, mertens_y,
                       primes, pi, sieve, events, sorter, max_m);
}

uint128_t D_hyperbola(uint128_t xp,
                      int64_t p,
                      int64_t b,
                      int64_t y,
                      int64_t low,
                      int64_t high,
                      int64_t phi,
                      int64_t mertens_y,
                      const Vector<int64_t>& primes,
                      const PiTable& pi,
                      Sieve& sieve,
                      Vector<SemiprimeEvent>& events,
                      SemiprimeEventSorter& sorter,
                      int64_t& max_m)
{
  return ::D_hyperbola(xp, p, b, y, low, high, phi, mertens_y,
                       primes, pi, sieve, events, sorter, max_m);
}

#endif

} // namespace
