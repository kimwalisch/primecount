///
/// @file  DBand.hpp
/// @brief Cancel prime leaves in the central regions of Gourdon's D
///        formula using a weighted Dirichlet hyperbola decomposition.
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

struct CentralEvent
{
  int64_t n;
  int64_t weight;
};

/// For p = primes[b] and p^2 <= t < p^3, the signed rough-number
/// prefix M(t) = sum_{n <= t, P-(n) > p} mu(n) satisfies
/// M(t) = phi(t, b - 1) - central_correction(t, p, b).
///
inline int64_t central_correction(int64_t t,
                                  int64_t p,
                                  int64_t b,
                                  const PiTable& pi)
{
  return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * b + 2;
}

/// Sum the correction over prime n in (lower, upper]. The last two
/// terms are inverted sums, avoiding two pi lookups per prime n.
///
template <typename T, typename Primes>
T central_prime_sum(T xp,
                    int64_t p,
                    int64_t b,
                    int64_t lower,
                    int64_t upper,
                    const Primes& primes,
                    const PiTable& pi)
{
  int64_t begin = pi[lower];
  int64_t end = pi[upper];
  int64_t count = end - begin;
  if (!count)
    return 0;

  T sum = 0;
  int64_t i = begin + 1;
  for (; i + 3 <= end; i += 4)
  {
    int64_t t0 = fast_div64(xp, primes[i]);
    int64_t t1 = fast_div64(xp, primes[i + 1]);
    int64_t t2 = fast_div64(xp, primes[i + 2]);
    int64_t t3 = fast_div64(xp, primes[i + 3]);
    sum += pi[t0] + pi[t1] + pi[t2] + pi[t3];
  }
  for (; i <= end; i++)
    sum += pi[fast_div64(xp, primes[i])];
  sum *= 2;

  T k = xp / p;
  int64_t r0 = fast_div64(k, upper);
  int64_t r1 = fast_div64(k, lower + 1);
  int64_t s0 = isqrt(fast_div64(xp, upper));
  int64_t s1 = isqrt(fast_div64(xp, lower + 1));
  sum += T(count) * (pi[r0] + pi[s0] - 4 * b + 2);
  for (int64_t j = pi[r0] + 1, last = pi[r1]; j <= last; j++)
    sum += pi[fast_div64(k, primes[j])] - begin;
  for (int64_t j = pi[s0] + 1, last = pi[s1]; j <= last; j++)
  {
    int64_t r = primes[j];
    sum += pi[fast_div64(xp, r * r)] - begin;
  }
  return sum;
}

/// Let X = floor(x/p), A = floor(X/U), H = floor(sqrt(X)) and
/// Phi(t) = phi(t, b - 1), with p^2 <= A <= H < U < p^3 and U <= y.
/// The contribution from A < m <= U is
///   M(H) Phi(H) - M(U) Phi(A)
///   + sum_{A < n <= H} a(n) E(X/n)
///   - sum_{A < n <= H} (a(n) + c(n)) Phi(X/n),
/// where a(n) = [P-(n) >= p], c(n) = mu(n)[P-(n) > p], E = Phi - M.
/// Prime n > p have a(n) + c(n) = 0. The remaining terms have n = q*r,
/// p <= q <= r, with weight 1 if q = p or q = r, and weight 2 otherwise.
/// Events are processed in decreasing n so sieve counts increase.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T central_events(T xp,
                               int64_t p,
                               int64_t b,
                               int64_t lower,
                               int64_t middle,
                               int64_t low,
                               int64_t high,
                               int64_t phi,
                               const Primes& primes,
                               const PiTable& pi,
                               const Count& count,
                               Vector<CentralEvent>& events,
                               SemiprimeEventBitmap& sorter)
{
  T sum = 0;
  if (low <= middle && middle < high)
  {
    int64_t value = phi + count(middle - low);
    sum += T(value) * (value - central_correction(middle, p, b, pi));
  }

  int64_t lower_n = min(fast_div(xp, high), middle);
  lower_n = max(lower_n, lower);
  int64_t upper_n = min(fast_div(xp, max(low, int64_t(1))), middle);
  if (lower_n < upper_n)
  {
    sum += central_prime_sum(xp, p, b, lower_n, upper_n, primes, pi);

    events.clear();
    int64_t last_q = pi[isqrt(upper_n)];
    for (int64_t i = b; i <= last_q; i++)
    {
      int64_t q = primes[i];
      int64_t first_r = max(i, pi[lower_n / q] + 1);
      int64_t last_r = pi[upper_n / q];
      for (int64_t j = last_r; j >= first_r; j--)
        events.push_back(CentralEvent{q * primes[j], (i == b || i == j) ? 1 : 2});
    }
    sorter.sort(events, p, lower_n, upper_n);

    for (const CentralEvent& event : events)
    {
      int64_t t = fast_div64(xp, event.n);
      int64_t value = phi + count(t - low);
      sum += central_correction(t, p, b, pi);
      sum -= T(event.weight) * value;
    }
  }
  return sum;
}

template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T central_y_segment(T xp,
                                  int64_t p,
                                  int64_t b,
                                  int64_t lower,
                                  int64_t middle,
                                  int64_t mertens_y,
                                  int64_t low,
                                  int64_t high,
                                  int64_t phi,
                                  const Primes& primes,
                                  const PiTable& pi,
                                  const Count& count,
                                  Vector<CentralEvent>& events,
                                  SemiprimeEventBitmap& sorter)
{
  T sum = 0;
  if (low <= lower && lower < high)
  {
    int64_t value = phi + count(lower - low);
    sum -= T(value) * mertens_y;
  }
  sum += central_events(xp, p, b, lower, middle, low, high,
                        phi, primes, pi, count, events, sorter);
  return sum;
}

/// Rows with x/p^3 < y and sqrt(x/p) < y < p^3. In these rows the
/// band (x/(p*y), y] of D contains all hard leaves above x/(p*y) and
/// all easy leaves of the row, hence the C1 formula skips these rows.
///
template <typename T>
ALWAYS_INLINE bool is_extended_row(T xp, int64_t p, int64_t y)
{
  int64_t upper = min(fast_div(xp, p * p), y + 1);
  int64_t middle = isqrt(xp);
  return upper < y && p * p > y / p && middle < y;
}

/// Extended rows: the band terms whose positions are < p^2 only
/// need PiTable lookups. The segmented sieve never processes the
/// row p below p^2, hence these terms are computed once per row.
///
template <typename T, typename Primes>
T central_cheap_terms(T xp,
                      int64_t p,
                      int64_t b,
                      int64_t y,
                      int64_t mertens_y,
                      const Primes& primes,
                      const PiTable& pi)
{
  if (!is_extended_row(xp, p, y))
    return 0;

  int64_t square = p * p;
  int64_t upper = fast_div64(xp, square);
  int64_t lower = fast_div64(xp, y);
  int64_t middle = isqrt(xp);

  // lower = X/y < p^2, phi(lower, b - 1) counts
  // 1 and the primes in [p, lower].
  int64_t phi_lower = pi[lower] - b + 2;
  T sum = 0;
  sum -= T(mertens_y) * T(phi_lower);

  if (middle < square)
  {
    int64_t phi_middle = pi[middle] - b + 2;
    int64_t mu_middle = 1 - (pi[middle] - b);
    sum += T(mu_middle) * T(phi_middle);

    // Primes n in (X/p^2, H] have X/n < p^2:
    // E(X/n) = phi(X/n, b - 1) - M(X/n) = 2 * pi[X/n] - 2b + 1.
    int64_t begin = pi[max(upper, lower)];
    int64_t end = pi[middle];
    T pi_sum = 0;
    for (int64_t i = begin + 1; i <= end; i++)
      pi_sum += pi[fast_div64(xp, primes[i])];
    sum += pi_sum * 2;
    sum -= T(end - begin) * T(2 * b - 1);
  }

  return sum;
}

/// Extended rows, terms of the current segment [low, high[:
/// the band terms with positions >= p^2 and the remaining
/// hard leaves m <= X/y, which are prime since X/y < p^2.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE T central_extended_segment(T xp,
                                         int64_t p,
                                         int64_t b,
                                         int64_t y,
                                         int64_t low,
                                         int64_t high,
                                         int64_t phi,
                                         const Primes& primes,
                                         const PiTable& pi,
                                         const Count& count,
                                         Vector<CentralEvent>& events,
                                         SemiprimeEventBitmap& sorter)
{
  int64_t square = p * p;
  int64_t upper = fast_div64(xp, square);
  int64_t lower = fast_div64(xp, y);
  int64_t middle = isqrt(xp);
  T sum = 0;

  if (middle >= square)
    sum += central_events(xp, p, b, lower, middle, low, high,
                          phi, primes, pi, count, events, sorter);
  else
  {
    // n <= X/p^2 < p^2: all n are prime, no events
    int64_t lower_n = max(fast_div(xp, high), lower);
    int64_t upper_n = min(fast_div(xp, max(low, int64_t(1))), upper);
    if (lower_n < upper_n)
      sum += central_prime_sum(xp, p, b, lower_n, upper_n, primes, pi);
  }

  // The remaining hard leaves m <= X/y are prime
  int64_t lower_m = min(fast_div(xp, high), lower);
  lower_m = max3(lower_m, y / p, p);
  int64_t upper_m = min(fast_div(xp, max(low, int64_t(1))), lower);
  if (lower_m < upper_m)
    for (int64_t i = pi[upper_m]; primes[i] > lower_m; i--)
    {
      int64_t t = fast_div64(xp, primes[i]);
      sum += phi + count(t - low);
    }
  return sum;
}

/// Process the central region of the current row p = primes[b].
/// Returns true if all leaves of this row have been processed.
/// Otherwise max_m may be lowered and the caller processes the
/// remaining leaves m <= max_m using the factor table.
///
template <typename T, typename Primes, typename Count>
ALWAYS_INLINE bool central_row(T xp,
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
                               Vector<CentralEvent>& events,
                               SemiprimeEventBitmap& sorter,
                               T& sum,
                               int64_t& max_m)
{
  int64_t upper = min(fast_div(xp, p * p), y + 1);

  // Rows with x/p^3 = y have no easy leaves,
  // they use the y-region with x/(p*y) = p^2.
  if (upper < y)
  {
    if (is_extended_row(xp, p, y))
    {
      sum += central_extended_segment(xp, p, b, y, low, high, phi,
                                      primes, pi, count, events, sorter);
      return true;
    }
  }
  else if (p * p > y / p)
  {
    int64_t middle = isqrt(xp);
    if (middle < y)
    {
      int64_t lower = fast_div64(xp, y);
      sum += central_y_segment(xp, p, b, lower, middle,
                               mertens_y[b], low, high, phi,
                               primes, pi, count, events, sorter);
      max_m = min(max_m, lower);
    }
  }

  return false;
}

} // namespace

#endif
