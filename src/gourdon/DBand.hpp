///
/// @file  DBand.hpp
/// @brief Cancel the prime leaves in the band (x/(p*y), y] of the rows
///        of Gourdon's D formula using a weighted Dirichlet hyperbola
///        rearrangement.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_BAND_HPP
#define D_BAND_HPP

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

  T sum = 0;
  int64_t i = pi_min_n + 1;
  for (; i + 3 <= pi_max_n; i += 4)
  {
    int64_t t0 = fast_div64(xp, primes[i]);
    int64_t t1 = fast_div64(xp, primes[i + 1]);
    int64_t t2 = fast_div64(xp, primes[i + 2]);
    int64_t t3 = fast_div64(xp, primes[i + 3]);
    sum += pi[t0] + pi[t1] + pi[t2] + pi[t3];
  }
  for (; i <= pi_max_n; i++)
    sum += pi[fast_div64(xp, primes[i])];
  sum *= 2;

  T xpp = xp / p;
  int64_t min_r = fast_div64(xpp, max_n);
  int64_t max_r = fast_div64(xpp, min_n + 1);
  int64_t min_r2 = isqrt(fast_div64(xp, max_n));
  int64_t max_r2 = isqrt(fast_div64(xp, min_n + 1));
  sum += T(pi_max_n - pi_min_n) * (pi[min_r] + pi[min_r2] - 4 * b + 2);
  for (int64_t j = pi[min_r] + 1, last = pi[max_r]; j <= last; j++)
    sum += pi[fast_div64(xpp, primes[j])] - pi_min_n;
  for (int64_t j = pi[min_r2] + 1, last = pi[max_r2]; j <= last; j++)
  {
    int64_t r = primes[j];
    sum += pi[fast_div64(xp, r * r)] - pi_min_n;
  }
  return sum;
}

/// Let X = floor(x/p), L = lower = floor(X/U), H = sqrt_xp =
/// floor(sqrt(X)) and Phi(t) = phi(t, b - 1), with p <= L <= H < U,
/// p^2 <= H and U <= y < p^3. The contribution from L < m <= U is
///   M(H) Phi(H) - M(U) Phi(L)
///   + sum_{L < n <= H} a(n) E(X/n)
///   - sum_{L < n <= H} (a(n) + c(n)) Phi(X/n),
/// where a(n) = [lpf(n) >= p], c(n) = mu(n)[lpf(n) > p], E = Phi - M.
/// This function computes all terms except -M(U) Phi(L) whose
/// positions are inside the current segment [low, high[.
/// Prime n > p have a(n) + c(n) = 0. The remaining terms have n = q*r,
/// p <= q <= r, with weight 1 if q = p or q = r, and weight 2 otherwise.
/// Events are processed in decreasing n so sieve counts increase.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T band_terms(T xp,
                           int64_t p,
                           int64_t b,
                           int64_t lower,
                           int64_t sqrt_xp,
                           int64_t low,
                           int64_t high,
                           int64_t phi,
                           const Primes& primes,
                           const PiTable& pi,
                           const Count& count,
                           Vector<SemiprimeEvent>& events,
                           SemiprimeEventSorter& sorter)
{
  T sum = 0;
  if (low <= sqrt_xp && sqrt_xp < high)
  {
    int64_t phi_sqrt_xp = phi + count(sqrt_xp - low);
    sum += T(phi_sqrt_xp) * (phi_sqrt_xp - phi_minus_mertens(sqrt_xp, p, b, pi));
  }

  int64_t min_n = min(fast_div(xp, high), sqrt_xp);
  min_n = max(min_n, lower);
  int64_t max_n = min(fast_div(xp, max(low, int64_t(1))), sqrt_xp);
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

/// Rows with x/p^3 >= y: all leaves of the band (x/(p*y), y]
/// are hard leaves, terms of the current segment [low, high[.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T hard_band_segment(T xp,
                                  int64_t p,
                                  int64_t b,
                                  int64_t lower,
                                  int64_t sqrt_xp,
                                  int64_t mertens_y,
                                  int64_t low,
                                  int64_t high,
                                  int64_t phi,
                                  const Primes& primes,
                                  const PiTable& pi,
                                  const Count& count,
                                  Vector<SemiprimeEvent>& events,
                                  SemiprimeEventSorter& sorter)
{
  T sum = 0;
  if (low <= lower && lower < high)
  {
    int64_t phi_lower = phi + count(lower - low);
    sum -= T(phi_lower) * mertens_y;
  }
  sum += band_terms(xp, p, b, lower, sqrt_xp, low, high,
                    phi, primes, pi, count, events, sorter);
  return sum;
}

/// Rows with x/p^3 < y and sqrt(x/p) < y < p^3. In these rows the
/// band (x/(p*y), y] of D contains all hard leaves above x/(p*y) and
/// all easy leaves of the row, hence the C1 formula skips these rows.
///
template <typename T>
ALWAYS_INLINE bool band_has_easy_leaves(T xp, int64_t p, int64_t y)
{
  int64_t x_div_p3 = min(fast_div(xp, p * p), y + 1);
  int64_t sqrt_xp = isqrt(xp);
  return x_div_p3 < y && p * p > y / p && sqrt_xp < y;
}

/// Rows whose band has easy leaves: the band terms whose positions
/// are < p^2 only need PiTable lookups. D_thread() skips the row p
/// in threads whose sieve limit is < p^2, hence these terms are
/// computed once per row and mixed_band_segment() only computes the
/// band terms whose positions are >= p^2.
///
template <typename T, typename Primes>
T easy_band_terms(T xp,
                  int64_t p,
                  int64_t b,
                  int64_t y,
                  int64_t mertens_y,
                  const Primes& primes,
                  const PiTable& pi)
{
  if (!band_has_easy_leaves(xp, p, y))
    return 0;

  int64_t square = p * p;
  int64_t x_div_p3 = fast_div64(xp, square);
  int64_t lower = fast_div64(xp, y);
  int64_t sqrt_xp = isqrt(xp);

  // lower = X/y < p^2, phi(lower, b - 1) counts
  // 1 and the primes in [p, lower].
  int64_t phi_lower = pi[lower] - b + 2;
  T sum = 0;
  sum -= T(mertens_y) * T(phi_lower);

  if (sqrt_xp < square)
  {
    int64_t phi_sqrt_xp = pi[sqrt_xp] - b + 2;
    int64_t mertens_sqrt_xp = 1 - (pi[sqrt_xp] - b);
    sum += T(mertens_sqrt_xp) * T(phi_sqrt_xp);

    // Primes n in (X/p^2, H] have X/n < p^2:
    // E(X/n) = phi(X/n, b - 1) - M(X/n) = 2 * pi[X/n] - 2b + 1.
    int64_t pi_min_n = pi[max(x_div_p3, lower)];
    int64_t pi_max_n = pi[sqrt_xp];
    T pi_sum = 0;
    for (int64_t i = pi_min_n + 1; i <= pi_max_n; i++)
      pi_sum += pi[fast_div64(xp, primes[i])];
    sum += pi_sum * 2;
    sum -= T(pi_max_n - pi_min_n) * T(2 * b - 1);
  }

  return sum;
}

/// Rows whose band has easy leaves, terms of the current segment
/// [low, high[: the band terms with positions >= p^2 and the
/// remaining hard leaves m <= X/y, which are prime since X/y < p^2.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T mixed_band_segment(T xp,
                                   int64_t p,
                                   int64_t b,
                                   int64_t y,
                                   int64_t low,
                                   int64_t high,
                                   int64_t phi,
                                   const Primes& primes,
                                   const PiTable& pi,
                                   const Count& count,
                                   Vector<SemiprimeEvent>& events,
                                   SemiprimeEventSorter& sorter)
{
  int64_t square = p * p;
  int64_t x_div_p3 = fast_div64(xp, square);
  int64_t lower = fast_div64(xp, y);
  int64_t sqrt_xp = isqrt(xp);
  T sum = 0;

  if (sqrt_xp >= square)
    sum += band_terms(xp, p, b, lower, sqrt_xp, low, high,
                      phi, primes, pi, count, events, sorter);
  else
  {
    // n <= X/p^2 < p^2: all n are prime, no events
    int64_t min_n = max(fast_div(xp, high), lower);
    int64_t max_n = min(fast_div(xp, max(low, int64_t(1))), x_div_p3);
    if (min_n < max_n)
      sum += sum_phi_minus_mertens(xp, p, b, min_n, max_n, primes, pi);
  }

  // The remaining hard leaves m <= X/y are prime
  int64_t min_m = min(fast_div(xp, high), lower);
  min_m = max3(min_m, y / p, p);
  int64_t max_m = min(fast_div(xp, max(low, int64_t(1))), lower);
  if (min_m < max_m)
    for (int64_t i = pi[max_m]; primes[i] > min_m; i--)
    {
      int64_t xpq = fast_div64(xp, primes[i]);
      sum += phi + count(xpq - low);
    }
  return sum;
}

/// Process the band of the current row p = primes[b].
/// Returns true if all leaves of this row have been processed.
/// Otherwise max_m may be lowered and the caller processes the
/// remaining leaves m <= max_m using the factor table.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE bool band_row(T xp,
                            int64_t p,
                            int64_t b,
                            int64_t y,
                            int64_t low,
                            int64_t high,
                            int64_t phi,
                            const Primes& primes,
                            const PiTable& pi,
                            const Vector<int64_t>& mertens_y,
                            const Count& count,
                            Vector<SemiprimeEvent>& events,
                            SemiprimeEventSorter& sorter,
                            T& sum,
                            int64_t& max_m)
{
  int64_t x_div_p3 = min(fast_div(xp, p * p), y + 1);

  // Rows with x/p^3 = y have no easy leaves,
  // their band starts at x/(p*y) = p^2.
  if (x_div_p3 < y)
  {
    if (band_has_easy_leaves(xp, p, y))
    {
      sum += mixed_band_segment(xp, p, b, y, low, high, phi,
                                primes, pi, count, events, sorter);
      return true;
    }
  }
  else if (p * p > y / p)
  {
    int64_t sqrt_xp = isqrt(xp);
    if (sqrt_xp < y)
    {
      int64_t lower = fast_div64(xp, y);
      sum += hard_band_segment(xp, p, b, lower, sqrt_xp,
                               mertens_y[b], low, high, phi,
                               primes, pi, count, events, sorter);
      max_m = min(max_m, lower);
    }
  }

  return false;
}

} // namespace

#endif
