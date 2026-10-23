///
/// @file  DHyperbola.hpp
/// @brief Compute the special leaves x/(p*m) with x/(p*y) < m <= y of
///        Gourdon's D formula using a weighted Dirichlet hyperbola
///        rearrangement, in which the prime leaves cancel.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_HYPERBOLA_HPP
#define D_HYPERBOLA_HPP

#include "SemiprimeEventSorter.hpp"

#include <fast_div.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {

struct SemiprimeEvent
{
  int64_t n;
  int64_t weight;
};

/// For p = primes[b] and p^2 <= t < p^3, the Moebius prefix sum
/// over the p-rough integers M(t) = sum_{n <= t, lpf(n) > p} mu(n)
/// satisfies M(t) = phi(t, b - 1) - phi_minus_mertens(t, p, b, pi).
///
inline int64_t phi_minus_mertens(int64_t t,
                                 int64_t p,
                                 int64_t b,
                                 const PiTable& pi)
{
  return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * b + 2;
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

/// Sum phi_minus_mertens(x / (p * q)) over the primes q in
/// (min_n, max_n]. Gourdon's inversion equality turns the sums of
/// pi(x / (p^2 * q)) and pi(sqrt(x / (p * q))) into sums over the
/// reflected primes r, avoiding two pi lookups per prime q.
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
  int64_t min_r2 = isqrt(fast_div64(xp, max_n));
  int64_t max_r2 = isqrt(fast_div64(xp, min_n + 1));

  T sum = sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
  sum += sum_pi(xpp, pi[min_r] + 1, pi[max_r], primes, pi);
  sum += T(pi_max_n) * pi[min_r] - T(pi_min_n) * pi[max_r];
  sum += T(pi_max_n) * pi[min_r2] - T(pi_min_n) * pi[max_r2];
  sum -= T(pi_max_n - pi_min_n) * (4 * b - 2);

  for (int64_t j = pi[min_r2] + 1, last = pi[max_r2]; j <= last; j++)
  {
    int64_t r = primes[j];
    sum += pi[fast_div64(xp, r * r)];
  }

  return sum;
}

/// Primes (x/y)^(1/3) < p <= sqrt(y) have x < p^3 * y <= p * y^2
/// and y^2 <= x, hence sqrt(x/p) < y < p^3. So hyperbola_sum()
/// computes their leaves x/(p*y) < m <= y, which also include their
/// easy leaves. The hyperbola terms whose positions are < p^2 only
/// need PiTable lookups. D_thread() skips the prime p in threads
/// whose sieve limit is < p^2, hence these terms are computed once
/// per prime and hyperbola_sum() only computes the hyperbola terms
/// whose positions are >= p^2.
///
template <typename T, typename Primes>
T easy_hyperbola_sum(T xp,
                     int64_t p,
                     int64_t b,
                     int64_t y,
                     int64_t mertens_y,
                     const Primes& primes,
                     const PiTable& pi)
{
  int64_t square = p * p;
  int64_t x_div_p3 = fast_div64(xp, square);
  int64_t xpy = fast_div64(xp, y);
  int64_t sqrt_xp = isqrt(xp);

  // xpy = X/y < p^2, phi(xpy, b - 1) counts
  // 1 and the primes in [p, xpy].
  int64_t phi_xpy = pi[xpy] - b + 2;
  T sum = 0;
  sum -= T(mertens_y) * T(phi_xpy);

  if (sqrt_xp < square)
  {
    int64_t phi_sqrt_xp = pi[sqrt_xp] - b + 2;
    int64_t mertens_sqrt_xp = 1 - (pi[sqrt_xp] - b);
    sum += T(mertens_sqrt_xp) * T(phi_sqrt_xp);

    // Primes n in (X/p^2, H] have X/n < p^2 (X/p^2 >= L as p^2 <= y):
    // E(X/n) = phi(X/n, b - 1) - M(X/n) = 2 * pi[X/n] - 2b + 1.
    int64_t pi_min_n = pi[x_div_p3];
    int64_t pi_max_n = pi[sqrt_xp];
    sum += sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
    sum -= T(pi_max_n - pi_min_n) * T(2 * b - 1);
  }

  return sum;
}

/// Compute the special leaves x/(p*y) < m <= y of the prime
/// p = primes[b] in the current segment [low, high[ using the
/// hyperbola method, if sqrt(x/p) < y < p^3. Then max_m is lowered
/// to x/(p*y) and the caller processes the remaining leaves using
/// the factor table.
///
/// Let X = floor(x/p), L = xpy = floor(X/y), H = sqrt_xp =
/// floor(sqrt(X)) and Phi(t) = phi(t, b - 1). The leaves L < m <= y
/// contribute
///   M(H) Phi(H) - M(y) Phi(L)
///   + sum_{L < n <= H} a(n) E(X/n)
///   - sum_{L < n <= H} (a(n) + c(n)) Phi(X/n),
/// where a(n) = [lpf(n) >= p], c(n) = mu(n)[lpf(n) > p], E = Phi - M.
/// Here only the terms whose positions are >= p^2 are computed, see
/// easy_hyperbola_sum(). Prime n > p have a(n) + c(n) = 0. The
/// remaining terms have n = q*r, p <= q <= r, with weight 1 if q = p
/// or q = r, and weight 2 otherwise. Events are processed in
/// decreasing n so sieve counts increase.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T hyperbola_sum(T xp,
                              int64_t p,
                              int64_t b,
                              int64_t y,
                              int64_t low,
                              int64_t high,
                              int64_t phi,
                              int64_t mertens_y,
                              const Primes& primes,
                              const PiTable& pi,
                              const Count& count,
                              Vector<SemiprimeEvent>& events,
                              SemiprimeEventSorter& sorter,
                              int64_t& max_m)
{
  int64_t square = p * p;
  int64_t sqrt_xp = isqrt(xp);
  if (sqrt_xp >= y || square <= y / p)
    return 0;

  T sum = 0;
  int64_t xpy = fast_div64(xp, y);
  max_m = min(max_m, xpy);
  int64_t query_low = max(low, square);

  if (query_low <= xpy && xpy < high)
  {
    int64_t phi_xpy = phi + count(xpy - low);
    sum -= T(phi_xpy) * mertens_y;
  }

  if (query_low <= sqrt_xp && sqrt_xp < high)
  {
    int64_t phi_sqrt_xp = phi + count(sqrt_xp - low);
    sum += T(phi_sqrt_xp) * (phi_sqrt_xp - phi_minus_mertens(sqrt_xp, p, b, pi));
  }

  // n in (L, H] with position X/n inside [max(low, p^2), high[
  int64_t min_n = min(fast_div(xp, high), sqrt_xp);
  min_n = max(min_n, xpy);
  int64_t max_n = min(fast_div(xp, query_low), sqrt_xp);

  if (min_n < max_n)
  {
    sum += sum_phi_minus_mertens(xp, p, b, min_n, max_n, primes, pi);

    events.clear();
    int64_t max_i = pi[isqrt(max_n)];
    for (int64_t i = b; i <= max_i; i++)
    {
      int64_t q = primes[i];
      int64_t min_j = max(i, pi[min_n / q] + 1);
      int64_t max_j = pi[max_n / q];
      for (int64_t j = max_j; j >= min_j; j--)
        events.push_back(SemiprimeEvent{q * primes[j], (i == b || i == j) ? 1 : 2});
    }
    sorter.sort(events, p, min_n, max_n);

    for (const SemiprimeEvent& event : events)
    {
      int64_t xpn = fast_div64(xp, event.n);
      int64_t phi_xpn = phi + count(xpn - low);
      sum += phi_minus_mertens(xpn, p, b, pi);
      sum -= T(event.weight) * phi_xpn;
    }
  }

  return sum;
}

} // namespace

#endif
