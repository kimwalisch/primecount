# The hyperbola method for the hard special leaves of Gourdon's D formula

Specification of the "full extension" and its implementation in primecount

- Repository: `kimwalisch/primecount`, branch `claude/d-extended-band`,
  head `22f2be25` (2026-10-23), after the C1 removal (Section 10.3), the
  reduction of AC's primes vector (Section 4.5) and the switch to the
  hyperbola terminology (Section 2).
- Base: PR #146 by Chopper Lin (head `cc9e282`), which is based on
  master `bcbf2f6`.
- Companion document: `Hard_Special_Leaves_Cancellation_2.tex` (Kim
  Walisch, work in progress). It derives the method of PR #146. This
  document describes the extension of that method and the current code.

---

## 0. How to use this document

The tex paper explains the mathematics of PR #146:
- the functions Phi_p, M_p and E_p (paper Section 3);
- the hyperbola identity (Section 4, called the band identity there);
- the cancellation of the prime leaves (Section 5);
- the two regions of PR #146, the *natural region* and the *y-region*
  (Sections 6 and 7).

The code on this branch goes further:

1. Every prime with `sqrt(x/p) < y < p^3` (and `p <= sqrt(y)`) computes
   its special leaves `x/(p*y) < m <= y` with the hyperbola method. If
   `x/p^3 < y`, these also include the easy leaves of `p`, which
   Gourdon's formula C (part C1) used to compute. The natural region of
   PR #146 no longer exists in the code.
2. The code was renamed and restructured after PR #146. It uses the
   established terms hyperbola method, special leaves and prime `p`
   (Section 2), whereas the PR and the paper use region, band and row.

Where the paper and this document disagree about regions or names,
this document describes the code.

Your tasks:

- **A. Full code review.** Verify every formula and bound against the
  mathematics; Section 9 has a checklist and Sections 3–5 the proofs
  to verify.
- **B. Review the C1 removal.** C1 was empty for every valid `y` and
  has been deleted without changing `y` or `alpha` (Section 10).
- **C. Further refinements.** Section 11 lists candidates.

Read `AGENTS.md` in the repository root first. It covers coding
conventions, required paper reading, ChangeLog rules and required
verification. Also run a **Debug build** (assertions enabled), as CI does
(Section 6). A bug in this extension was found by CI only because Release
builds compile `ASSERT()` out.

---

## 1. Summary

primecount computes `pi(x) = A - B + C + D + Phi0 + Sigma` (Gourdon) with
`z = y`. Formula D sums the hard special leaves:

```
D = - sum_{k < b <= pi(x^(1/4))} sum_m mu(m) phi(x/(p_b m), b - 1)
```

The sum runs over the square-free `m` with `lpf(m) > p`, `m <= y < p*m`
and `m <= x/p^3`. The C formula sums the easy leaves (`m > x/p^3`) with
`phi(x/(pm), b-1) = pi(x/(pm)) - b + 2`. AC computes them in two parts:
C1 for `(x/y)^(1/3) < p <= sqrt(y)` and C2 for `sqrt(y) < p <= x_star`.

For a prime `p = p_b`, let `X = floor(x/p)`, `H = floor(sqrt(X))` and
`L = floor(X/y)`. A weighted Dirichlet hyperbola rearrangement (paper
Section 4, with `U = y`) rewrites the sum over the special leaves
`L < m <= y` of `p` as:

- two boundary terms;
- a cheap term `E_p(X/n)` for every `p`-rough `n` in `(L, H]`;
- one sieve query per semiprime `n` in `(L, H]`, the "events".

We call these the *hyperbola terms* of `p`. The sieve queries of all
prime `n` cancel. The position of every hyperbola term (`L`, `H` or
`X/n`) is `< y`.

A prime `p <= sqrt(y)` uses the hyperbola method if and only if
`sqrt(x/p) < y < p^3`, i.e. `p > max(y^(1/3), x/y^2)`. There are two
cases:

- **Hard leaves only** (`x/p^3 >= y`, equivalently `L >= p^2`): every
  leaf `L < m <= y` is a hard leaf of D. The remaining leaves `m <= L`
  are computed as before, using the factor table.
- **Also easy leaves** (`x/p^3 < y`, equivalently `L < p^2`,
  equivalently `p > (x/y)^(1/3)`): the leaves `L < m <= y` also include
  the easy leaves `m in (x/p^3, y]` of `p`, which are exactly C1's
  leaves of `p`. D computes them with the hyperbola method, and AC has no
  C1 formula. The remaining leaves `m <= L < p^2` are primes and are
  processed by the D kernel's normal factor-table loop.

Hyperbola terms whose position is `< p^2` need only PiTable lookups (the
*easy hyperbola terms*). They are computed once per prime by
`easy_hyperbola_terms()`. All other hyperbola terms are computed per
sieve segment by `hyperbola_sum()`.

For every valid `y`, every prime of C1, `(x/y)^(1/3) < p <= sqrt(y)`,
uses the hyperbola method and has easy leaves. Hence C1 computed nothing,
and it has been deleted without changing `y` or `alpha` (Section 10).
Primes that use the hyperbola method, and C1 primes, exist only if
`y > x^(2/5)`. At default parameters, D uses the hyperbola method only
for `x >= ~4*10^11`.

Effect: the share of D's sieve queries that are removed grows by 3–7
percentage points (Section 7). For example, at `x = 10^20` it rises from
31.4% (PR #146) to 36.6%, and at `x = 10^24` from 47.5% to 54.1%. AC gets
faster (C1 disappears) and D is unchanged within noise.

---

## 2. Notation and mapping to the paper and the code

| Math (paper) | Code | Meaning |
|---|---|---|
| `p = p_b`, `b` | `p` / `prime`, `b` | sieving prime and its index, `b = pi(p)` |
| `X = floor(x/p)` | `xp` | |
| `H = floor(sqrt(X))` | `sqrt_xp` | |
| `L = floor(X/U)` with `U = y` | `xpy` | the leaves `L < m <= y` use the hyperbola method |
| `U` | `y` | upper end of these leaves, the same for every prime |
| `floor(X/p^2) = floor(x/p^3)` | `x_div_p3`, `fast_div(xp, square)` | largest hard leaf; `n <= X/p^2` ⟺ position `X/n >= p^2` |
| `K = floor(X/p)` (paper Section 8) | `xpp` | `x/p^2` |
| `Phi_p(t) = phi(t, b-1)` | `phi + count(t - low)` (sieve) or `pi[t] - b + 2` (`t < p^2`) | |
| `M_p(t) = sum_{n<=t, lpf(n)>p} mu(n)` | "mertens" | Möbius prefix sum over the `p`-rough integers |
| `M_p(y)` | `mertens_y[b]` | precomputed once per prime |
| `E_p(t) = Phi_p(t) - M_p(t)` | `phi_minus_mertens(t, p, b, pi)` | valid for `p^2 <= t < p^3` |
| `a_p(n) = [lpf(n) >= p]`, `c_p(n) = mu(n)[lpf(n) > p]` | (comments) | |
| hyperbola terms | `hyperbola_sum()`, `easy_hyperbola_terms()` | the rearranged terms of the leaves `L < m <= y` |
| semiprime events | `SemiprimeEvent{n, weight}`, `SemiprimeEventSorter` | `n = q*r`, `p <= q <= r`, weight `w(q,r)` |
| `alpha < q <= beta` (Gourdon's inversion) | `min_n < q <= max_n` | |
| reflected primes `r` | `(min_r, max_r]`, `(min_r2, max_r2]` | |
| `m_max` | `max_m` | largest leaf of `p` in the segment |
| natural region (PR #146) | — | replaced by the hyperbola method with `U = y` |
| y-region (PR #146), "hard band" (earlier commits) | hard leaves only | primes with `x/p^3 >= y` |
| "extended rows", "mixed band" (earlier commits) | also easy leaves | primes with `x/p^3 < y` |
| "band", "row" (PR #146, the paper, earlier commits) | hyperbola terms, prime `p` | |
| "cheap terms", "easy band terms" (earlier commits) | easy hyperbola terms | positions `< p^2` |

The names follow the literature: Dirichlet's hyperbola method, and the
special leaves and easy and hard leaves of Lagarias–Miller–Odlyzko and
Deléglise–Rivat. Like Gourdon, who writes D as a sum of terms `(p, m)`,
the code refers to the special leaves of a prime `p` and needs no word
for the set of leaves of `p`. "Band" and "row" are not used in the
literature.

Least-prime-factor notation: the paper writes `P^-(n)`, the code comments
write `lpf(n)`, with `lpf(1) = +infinity`.

---

## 3. Mathematics

### 3.1 The functions Phi_p, M_p and E_p

These are taken from paper Section 3, with two additions for arguments
below `p^2`.

For `p^2 <= t < p^3`, every `n <= t` with `lpf(n) >= p` has at most two
prime factors (counted with multiplicity). The cancellation table of the
paper gives:

```
E_p(t) = Phi_p(t) - M_p(t) = 2 pi(t) + pi(t/p) + pi(sqrt(t)) - 4b + 2      (p^2 <= t < p^3)
```

`phi_minus_mertens()` implements this formula. It fails for `t < p^2`
(the `q^2` count becomes negative) and for `t = p^3`.

For `p <= t < p^2`, the `p`-rough integers `<= t` are 1 and the primes in
`(p, t]` (any composite with `lpf > p` exceeds `p^2`). The integers
counted by `Phi_p` are 1 and the primes in `[p, t]`. Hence:

```
Phi_p(t) = pi(t) - b + 2
M_p(t)   = 1 - (pi(t) - b)
E_p(t)   = 2 pi(t) - 2b + 1                                                (p <= t < p^2)
```

`easy_hyperbola_terms()` uses these three formulas, inline.

### 3.2 The hyperbola identity with U = y

This follows paper Section 4. The derivation (Steps 1–5) holds for any
weights. Step 6 needs every `t = X/n` with `L < n <= H` to satisfy
`p <= t < p^3`, and uses `E_p` in the form valid for that `t` (Section
3.1).

For a prime `p` with `H < y < p^3`, let `L = floor(X/y)`. Then:

```
D_p(L, y) := - sum_{L < m <= y} c_p(m) Phi_p(X/m)
           =  M_p(H) Phi_p(H) - M_p(y) Phi_p(L)
            + sum_{L < n <= H} a_p(n) E_p(X/n)
            - sum_{L < n <= H} (a_p(n) + c_p(n)) Phi_p(X/n)
```

Facts used by the code (each needs verification):

- `L <= H < y` because `H < y` ⟺ `X < y^2` ⟹ `X/y < y`, and `X/y <= sqrt(X)`.
- `L >= p`. We have `p <= x^(1/4)` (D's primes) and `y <= sqrt(x)`. So
  `p^2 * y <= sqrt(x) * sqrt(x) = x`, i.e. `x/(p*y) >= p`. Hence `n = p`
  never lies in `(L, H]`, and every prime `n` in `(L, H]` has
  `a + c = 1 - 1 = 0`.
- Positions: `X/n` for `L < n <= H` lies in `[H, X/(L+1)]`, and
  `X/(L+1) < y`. Hence the position of every hyperbola term (`L`, `H`,
  `X/n`) is in `[L, y)`, which is `>= p` and `< y < p^3`. So `E_p` is
  defined at every position, with the formula chosen by whether the
  position is `< p^2`.
- `M_p(y)` uses `y < p^3`, and `y >= p^2` holds because `p <= sqrt(y)`.
  Hence `M_p(y) = phi(y, b-1) - E_p(y)` with the `p^2 <= t < p^3`
  formula (Section 3.9).
- Weights for `L < n <= H` with `a_p(n) = 1`. Here `n <= H < y < p^3`, so
  `n` has at most two prime factors `>= p`:
  - prime `q > p`: weight 0, an `E_p` term only;
  - `p^2` and `p*q` (`q > p`): weight 1;
  - `q^2` (`q > p`): weight 1;
  - `q*r` (`p < q < r`): weight 2.

  When `L < p^2 <= H` (easy leaves), `n = p^2` is an event with weight 1.
  The event generator handles it (`i == b && j == b`).

### 3.3 Which leaves the hyperbola method computes

Every square-free `p`-rough `m` in `(L, y]` is a special leaf of `p`
(`m > p`, `m <= y`, `p*m > y`, largest prime factor `<= y`):

- `m > L >= p`, so `m > p`.
- `L = floor(x/(p*y)) >= floor(y/p)`, because `y/p <= x/(p*y)` ⟺
  `y^2 <= x`. So `m >= L + 1 > y/p`, i.e. `p*m > y = z`.
- `m <= y` ⟹ the largest prime factor of `m` is `<= y`.

A leaf `m` is **hard** if `m <= x/p^3` and **easy** otherwise:

- Hard leaves only (`x/p^3 >= y`): every leaf `L < m <= y` is hard.
  Then `L = floor(X/y) >= p^2`, so all hyperbola terms have positions
  `>= p^2`.
- Also easy leaves (`x/p^3 < y`): the leaves `L < m <= y` are the hard
  leaves `(L, x/p^3]` and the easy leaves `(x/p^3, y]`. The easy leaves
  are exactly C1's leaves of `p`. C1's range is
  `max(x/p^3, y/p) < m <= y`, and `x/p^3 >= y/p` because
  `p <= sqrt(y) <= sqrt(x/y)`. C1 adds `-mu(m) (pi(X/m) - b + 2)`, which
  is the same value as `-c_p(m) Phi_p(X/m)` for an easy leaf. Moving them
  into D changes A+C and D by the same `delta` with opposite signs.
  The unit-test values were regenerated as `AC_old - delta` and
  `D_old + delta`, with `delta` computed independently by brute force.

If `x/p^3 < y`, the remaining leaves `m in (y/p, L]` of `p` are hard
leaves. Since `L < p^2`, they are primes. D's factor-table loop processes
them (Section 4.3).

### 3.4 Classification of the primes

All bounds below are integer-exact, and all comparisons use floors.

| Condition in `hyperbola_sum()` | Equivalent | Proof sketch |
|---|---|---|
| `sqrt_xp < y` | `p > floor(x/y^2) = floor(floor(x/y)/y)` | `floor(sqrt(N)) < y` ⟺ `N < y^2`; `floor(x/p) < y^2` ⟺ `x < p y^2` ⟺ `floor(x/y^2) < p` |
| `p*p > y/p` (integer `y/p`) | `p > floor(y^(1/3))`, i.e. `y < p^3` | `floor(y/p) < p^2` ⟺ `y/p < p^2` |
| `xpy < square`, i.e. `floor(X/y) < p^2` | `x/p^3 < y` ⟺ `p > floor((x/y)^(1/3)) = iroot<3>(xy)` | `X/y < p^2` ⟺ `x/y < p^3` ⟺ `floor(x/y) < p^3` |

So, for `p <= sqrt(y)`, with the cutoff
`F = min(max(iroot<3>(y), xy / y), isqrt(y))`:

- `p` uses the hyperbola method ⟺ `p > F`.
- Its leaves `L < m <= y` also include easy leaves ⟺
  `p > iroot<3>(xy)`.

The `isqrt(y)` cap does not change the first equivalence for
`p <= sqrt(y)`, but it makes `F <= iroot<3>(xy)` hold for every valid `y`
(Section 10.1). The code does not compute `F` (it was `band_min_prime()`
until `60db7d7d`): `hyperbola_sum()` tests each prime itself, and `D.cpp`
uses `iroot<3>(xy)` (Section 4.4). All three conditions are monotone in
`p`, so each kind of prime forms a contiguous range of `b`. These
equivalences were checked by brute force on 259,957 random and edge cases
(Appendix A.1). The edge cases include `x = p^4, p^5, p^6 (±1)` and
`y = p^2, p^3 (±1)`.

Prime ranges, in increasing `p`:

| Primes | D | AC |
|---|---|---|
| `k < p <= F` | no hyperbola method; hard leaves as in master | (not a C1 prime) |
| `F < p <= (x/y)^(1/3)` | hyperbola method, hard leaves only | (not a C1 prime) |
| `(x/y)^(1/3) < p <= sqrt(y)` | hyperbola method, also easy leaves | C1 prime (C1 removed) |
| `sqrt(y) < p <= x^(1/4)` | 2-prime loop, no hyperbola method | C2 primes |

The ranges are in this order because `F <= iroot<3>(xy)` for every
valid `y` (Section 10.1); some of them may be empty. So D computes the
easy leaves of every C1 prime with the hyperbola method, which is why C1
could be deleted.

Primes that use the hyperbola method, and C1 primes, exist only if
`y > x^(2/5)`. The former have `x < p*y^2` (`H < y`), the latter
`x < p^3*y` (`p > (x/y)^(1/3)`), and with `p <= sqrt(y)` both give
`x < y^(5/2)`. For `y <= x^(2/5)`, every prime `p <= sqrt(y)` is in the
first range: D does not use the hyperbola method, and C1's range is
empty, already in master.

### 3.5 Position bookkeeping: easy hyperbola terms and segment terms

The D kernels process the prime `p = primes[b]` only in threads whose
sieve limit is `>= p^2`:
`max_b = pi[min3(isqrt(x/low1), isqrt(limit), x14)]` in `D_thread_*`.
Therefore the hyperbola terms at positions `< p^2` (whose `Phi_p` values
come from PiTable anyway) are computed once per prime in `D_OpenMP`.
`hyperbola_sum()` computes only the terms with positions `>= p^2`.

| Term | Position | Hard leaves only (`L >= p^2`) | Easy leaves, `H >= p^2` | Easy leaves, `H < p^2` |
|---|---|---|---|---|
| `-M_p(y) Phi_p(L)` | `L` | `hyperbola_sum()` | `easy_hyperbola_terms()` | `easy_hyperbola_terms()` |
| `M_p(H) Phi_p(H)` | `H` | `hyperbola_sum()` | `hyperbola_sum()` | `easy_hyperbola_terms()` |
| prime `n`: `E_p(X/n)` | `X/n` | `hyperbola_sum()` | `hyperbola_sum()` (all `X/n >= H >= p^2`) | `n <= X/p^2`: `hyperbola_sum()`; `n > X/p^2`: `easy_hyperbola_terms()` |
| events `n = q r` | `X/n` | `hyperbola_sum()` | `hyperbola_sum()` | none (all `n < p^2` are prime) |

`hyperbola_sum()` restricts all positions to `[query_low, high)` with
`query_low = max(low, p^2)`: `n <= X/query_low` for the prime terms and
events, `query_low <= L` and `query_low <= H` for the two boundary terms.
`easy_hyperbola_terms()` computes the complement.

**Proof that every per-segment term is computed exactly once.** Take a
position `t >= p^2`, `t < y`, in segment `[low, high)` of a thread with
limit `>= high > t`. Check the kernel's conditions for `b`:

1. `b <= max_b`:
   - `p^2 <= t < limit`, so `p <= isqrt(limit)`;
   - `low <= t < y <= X/p` (because `p^2 <= y <= sqrt(x)`), so
     `p^2 <= x/low`;
   - `p <= sqrt(y) <= x^(1/4)`.
2. `b >= min_b = max(k, pi[min(xy/limit, x14)]) + 1`. The hyperbola term
   positions are `>= L`, so `limit > t >= L` ⟹ `limit >= L+1 > x/(p*y)` ⟹
   `p > floor(xy/limit)`.
3. `prime >= max_m` (goto next segment) loses nothing. It triggers if
   `floor(x/p^3) <= p`, i.e. `p^4 <= x < p^3 (p+1)`. Then
   `X < p^3 + p^2`, `H < p^2` and `X/p^2 = p <= L`, so no `n` in `(L, H]`
   has a position `>= p^2`. No later prime `q > p` has `q^4 <= x` either.
   It also triggers if `X/low <= p`, i.e. `low >= X/p >= y`, but then no
   hyperbola term lies in the segment.
4. The prime-only cut-off (`max_m <= p^2 && low >= y`) requires
   `low >= y`, and all hyperbola term positions are `< y`.

### 3.6 Easy hyperbola terms (`easy_hyperbola_terms()`)

Only for primes with easy leaves (`x/p^3 < y`):

```
- M_p(y) Phi_p(L)                                with Phi_p(L) = pi(L) - b + 2 (p <= L < p^2)
+ M_p(H) Phi_p(H)                                if H < p^2, with the formulas of Section 3.1
+ sum_{q prime, X/p^2 < q <= H} (2 pi(X/q) - 2b + 1)    if H < p^2
```

The prime sum starts at `X/p^2`, which is `>= L` because `p^2 <= y`.

### 3.7 Prime terms with positions >= p^2 (`sum_phi_minus_mertens()`)

For primes `q` in `(alpha, beta] = (min_n, max_n]`, with all `X/q` in
`[p^2, p^3)`:

```
sum E_p(X/q) = 2 sum pi(X/q) + sum pi(floor(K/q)) + sum pi(floor(sqrt(X/q))) - N (4b - 2)
```

Here `K = floor(X/p) = xpp` and `N = pi(beta) - pi(alpha)`. Gourdon's
inversion equality, applied with floors, gives:

```
sum_{alpha < q <= beta} pi(floor(Z/q)) = pi(beta) pi(floor(Z/beta)) - pi(alpha) pi(floor(Z/(alpha+1)))
                                        + sum_{floor(Z/beta) < r <= floor(Z/(alpha+1))} pi(floor(Z/r))
```

Proof: count the prime pairs `(q, r)` with `q*r <= Z` and `q` in
`(alpha, beta]`, grouped by `r`:

- `r <= Z/beta`: all `N` values of `q` count.
- `Z/beta < r <= Z/(alpha+1)`: `q` runs over `(alpha, Z/r]`.
- `r > Z/(alpha+1)`: no `q`.

Collect the terms. The same argument with `q*r^2 <= X` handles
`pi(floor(sqrt(X/q)))`, with reflected primes `r` in
`(isqrt(X/beta), isqrt(X/(alpha+1))]` and summands `pi(floor(X/r^2))`.
`sum_phi_minus_mertens()` implements exactly this:

```
sum  = 2 * sum_pi(xp,  pi[min_n] + 1, pi[max_n])
     +     sum_pi(xpp, pi[min_r] + 1, pi[max_r])
     + pi[max_n] * pi[min_r]  - pi[min_n] * pi[max_r]
     + pi[max_n] * pi[min_r2] - pi[min_n] * pi[max_r2]
     - (pi[max_n] - pi[min_n]) * (4b - 2)
     + sum_{min_r2 < r <= max_r2} pi[X / r^2]
```

Here `min_r = K/max_n`, `max_r = K/(min_n+1)`,
`min_r2 = isqrt(X/max_n)` and `max_r2 = isqrt(X/(min_n+1))`. This is the
form used by C2 in `AC_*.hpp` (`pi[q_lo] * pi_y - pi[q_hi] * ...`).

### 3.8 M_p(y)

`D_OpenMP` precomputes `phi(y, b-1)` for all `b <= pi_sqrty` with
`phi_vector()`. It then subtracts `E_p(y) = phi_minus_mertens(y, ...)`.
The result is valid only for primes with `p^2 <= y < p^3`, which
includes all primes that use the hyperbola method. Other entries are
never read.

---

## 4. Implementation

### 4.1 Files

| File | Role |
|---|---|
| `src/gourdon/DHyperbola.hpp` | all hyperbola mathematics (252 lines) |
| `src/gourdon/SemiprimeEventSorter.hpp` | sorts events by decreasing `n` (bitmap over the mod-240 wheel, falls back to `std::sort`) — from PR #146 |
| `src/gourdon/D.cpp` | `D_OpenMP()`: PiTable/primes up to `y`, `mertens_y`, unsigned reduction, easy-terms loop |
| `src/gourdon/D_default.hpp`, `D_avx2.hpp`, `D_avx512.hpp`, `D_arm_neon.hpp`, `D_arm_sve.hpp` | D kernels, call `hyperbola_sum()` |
| `src/gourdon/AC.cpp` | AC's PiTable up to `y`, primes up to `min(isqrt(x / isqrt(y)), y)` (Section 4.5) |
| `src/gourdon/AC_default.hpp`, `AC_libdivide.hpp`, `AC_arm_sve.hpp` | A and C2 formulas; C1 was removed (Section 4.5) |
| `src/PiTable.hpp`, `PiTable.cpp` | `prev_prime(x)`: the largest prime `<= x`, used by AC |
| `test/gourdon/DHyperbola.cpp` | dense reference test of `hyperbola_sum()`/`easy_hyperbola_terms()` per segment, query order, `max_m`; exhaustive check that every prime `(x/y)^(1/3) < p <= sqrt(y)` uses the hyperbola method with easy leaves (`x <= 30000`) |
| `test/gourdon/SemiprimeEventSorter.cpp` | sorter test |
| `test/gourdon/D.cpp` | regenerated values; `D_reference()` brute force incl. the easy leaves computed by D; boundary cases; 30 hyperbola cases labelled 0 = not used, 1 = also easy leaves, 2 = hard leaves only |
| `test/gourdon/AC.cpp` | regenerated values (`AC_old - delta`) |

### 4.2 `DHyperbola.hpp`, function by function (line numbers at `22f2be25`)

- `struct SemiprimeEvent { int64_t n; int64_t weight; }` (line 29).
- `phi_minus_mertens(t, p, b, pi)` (line 39): `E_p(t)` for
  `p^2 <= t < p^3`.
- `sum_pi(xp, i, last, primes, pi)` (line 49): `sum_{i <= j <= last}
  pi[xp / primes[j]]`. It is unrolled 4x for ILP, and the tail loop is
  marked `NO_UNROLL_LOOP` (as in AC's `sum_pi`). Note that AC's
  `sum_pi<T, MULTIPLIER>` in `AC_default.hpp` is a *different* function in
  an anonymous namespace. Explicit template arguments keep the calls
  unambiguous.
- `sum_phi_minus_mertens(xp, p, b, min_n, max_n, primes, pi)` (line 80):
  Section 3.7. Precondition: `X/max_n >= p^2` and `X/(min_n+1) < p^3`.
- `easy_hyperbola_terms(xp, p, b, y, mertens_y, primes, pi)` (line 124):
  Section 3.6. Precondition: `p` uses the hyperbola method and has easy
  leaves. There is no internal test: `D.cpp` calls it for the primes
  `(x/y)^(1/3) < p <= sqrt(y)`, which all satisfy this (Section 10.1).
  Its doc comment gives the proof.
- `hyperbola_sum(xp, p, b, y, low, high, phi, mertens_y, primes, pi,
  count, events, sorter, max_m)` (line 180). Contract:
  - It tests the prime itself and returns 0 unless `H < y < p^3`
    (Section 3.4).
  - `phi + count(t - low) == Phi_p(t)` for `t` in `[low, high)`, and
    `count` must be called with nondecreasing arguments.
  - Returns the hyperbola terms with positions in `[max(low, p^2), high)`.
  - If `p` uses the hyperbola method, sets `max_m = min(max_m, L)`.
  - Order of sieve queries: `L` (if `L >= p^2`, i.e. hard leaves only),
    `H` (if `H >= p^2`),
    then events in decreasing `n` (increasing position `X/n >= H`).
  - The `n` range is computed as
    ```
    min_n = max(min(X/high, H), L)               // clamp with min() BEFORE narrowing to int64
    max_n = min(X/query_low, H)                  // query_low = max(low, p^2)
    ```
    so `n` in `(min_n, max_n]` ⟺ `X/n` in `[max(low, p^2), high)` and
    `n` in `(L, H]`. The `min(fast_div(xp, high), sqrt_xp)` clamp also
    prevents a 128-bit → `int64_t` narrowing; `max(fast_div(...), xpy)`
    would return the 128-bit type.
  - Events: for `q = primes[i]`, `i = b .. pi[isqrt(max_n)]`, the indices
    `j` run from `max(i, pi[min_n/q] + 1)` to `pi[max_n/q]`.
    `weight = (i == b || i == j) ? 1 : 2`.

### 4.3 Kernel integration (identical in all five D kernels)

```cpp
int64_t min_m = max(xp_high, y / prime);
int64_t max_m = min(fast_div(xp, prime * prime), xp_low);

if (prime >= max_m)
  goto next_segment;

// A composite m with all prime factors > prime exceeds
// prime^2 and all hyperbola terms have positions < y.
if (max_m <= prime * prime && low >= y)
  break;

sum += hyperbola_sum(xp, prime, b, y, low, high, phi[b], mertens_y[b],
                     primes, pi, count_sieve, events, sorter, max_m);

// ... unchanged factor table loop over (min_m, max_m] ...
```

- **Prime-only cut-off** (PR #146, Lin's report Section 7). The leaves of
  `p` in this segment satisfy `m <= max_m`, and `max_m` (uncapped)
  decreases in `p`. Once `max_m <= p^2`, every remaining leaf of `p` and
  of all larger primes is prime, so the 2-prime loop takes over. The test
  sits *before* `hyperbola_sum()`, so it uses the uncapped `max_m`. It
  requires `low >= y`, because the 2-prime loop computes no hyperbola
  terms; with `low >= y` there are none in the segment. Without it, the
  break would also skip the hyperbola terms of the current prime. This
  matters only if `p` uses the hyperbola method, which requires
  `y > x^(2/5)` (Section 3.4).
- **Remaining leaves.** For every prime that uses the hyperbola method,
  the factor-table loop computes the leaves `m <= min(max_m, L)`. Their
  positions `X/m >= X/L >= y` come after the positions of all hyperbola
  terms, so counts stay nondecreasing.
- In the 2-prime loop, `min_m = max3(xp_high, y / prime, prime)`. The
  `y / prime` bound is needed because primes with `b <= pi_sqrty` can
  reach this loop through the cut-off.
- `count_sieve` is a lambda (`sieve.count`), or the functors
  `SieveCountAvx512` / `SieveCountArmSve`. These must be functors because
  `operator()` needs the SIMD target attribute (see `AGENTS.md` on
  `PRIMECOUNT_TARGET_ATTRIBUTE`).
- `events` and `sorter` are created once per `D_thread_*` call and reused
  across primes and segments.

### 4.4 `D.cpp` (`D_OpenMP`)

- `PiTable pi(y)` and primes up to `y` (master had reduced them to
  `min(y, xy/isqrt(y))`). The hyperbola method needs `pi(t)` for `t < y`,
  `pi(y)` for `M_p(y)`, and primes up to `H < y`.
- `mertens_y = phi_vector(y, pi_sqrty)`, followed by
  `mertens_y[b] -= phi_minus_mertens(y, primes[b], b, pi)` for
  `b <= pi_sqrty`. All readers (the kernels' `hyperbola_sum()` calls and
  the easy-terms loop) have `b <= pi_sqrty`, and `y^2 <= x` implies
  `pi_sqrty <= pi[x14]`.
- Accumulation in unsigned `UT` (worker kernels and reduction).
  `D_thread<T>()` implicitly converts the unsigned worker result to signed
  `T` before assigning it to `thread.sum`. The final value is converted
  implicitly by the return from `D_OpenMP()`. The boundary products are
  bounded by `X` (`|M(y) Phi(L)| <= y L <= X`, `|M(H) Phi(H)| <= H^2 <= X`).
  Partial sums over primes and segments are not bounded tightly enough to
  rule out signed overflow (UB) for `T = int64_t` near `2^63`, so modular
  unsigned arithmetic plus a final signed representative is used. The
  final `D` is in signed range.
- Easy hyperbola terms: in the same OpenMP region, after a thread runs
  out of D work units, it processes
  `for (b = min_b++; b <= pi_sqrty; b = min_b++)` (`RelaxedAtomic`, as in
  `S2_easy.cpp`), with `min_b = max(k, pi_root3_xy) + 1`. These are
  exactly the primes with easy leaves (Section 3.4). Until `60db7d7d`,
  `min_b` also included
  `pi[band_min_prime(xy, y)]`, which is `<= pi_root3_xy` (Section 10.1).

### 4.5 AC

AC computes only the A and C2 formulas. `1aeaa05f` removed the C1
formula, i.e. the primes `max(k, pi_root3_xy) < b <= pi_sqrty` of C:
these are exactly D's primes with easy leaves (Section 3.4), and D
computes their easy leaves with the hyperbola method. The C2 doc comments
in `AC_*.hpp` refer to D for this part of C.

Before the removal, `AC_OpenMP_*` ran the C1 loop over
`b = min_c1 .. max_c1` with `max_c1 = pi[band_min_prime(xy, y)]` and
`min_c1 > max(k, pi_root3_xy)`. Since
`band_min_prime(xy, y) <= iroot<3>(xy)` (Section 10.1),
`max_c1 <= pi_root3_xy < min_c1`, so the loop was never entered, for any
valid `y`.

AC's PiTable is still sized `y`: C2's clustered leaves `q <= y` need
`pi(y)` and other lookups up to `y`. Since `bfb74102`, AC's primes vector
only goes up to `max_prime = min(isqrt(x / isqrt(y)), y)`:

- The primes `p` of A and C2 are `> sqrt(y)`, because
  `x_star >= x^(1/4) >= sqrt(y)` and C2 starts after `pi_sqrty`.
- They only read the primes `<= sqrt(x/p)`: the second primes of their
  sparse leaves, which are also `<= y` (C2's reflected range lies inside
  its sparse range), and the primes `p <= x^(1/3)` themselves, i.e.
  `p <= sqrt(x/p)`. All these primes `p` are `<= max_prime`; this was also
  checked for 95 million `(x, y)` pairs, exhaustively for `x <= 300000`
  and randomly up to `10^31`.
- C2's largest clustered prime `max_clustered_global`, the largest prime
  `<= y`, used to be `primes[pi[y]]` and now comes from
  `PiTable::prev_prime(y)`.

At default `y` the vector has 2.7x fewer primes at `1e19` and 4.8x fewer
at `1e24`, where they then fit `uint32_t`. AC's peak memory at `1e19`
dropped from 75.8 MB to 37.3 MB (libdivide stores a 9-byte divider per
prime). pi(x)'s peak memory is set by D (71.9 MB at `1e19`), which needs
primes up to almost `y`, because its smallest primes that use the
hyperbola method have `sqrt(x/p)` close to `y`.

### 4.6 `SemiprimeEventSorter`

The sorter comes from PR #146, renamed. Events must have distinct `n` in
`(min_n, max_n]`. Two bitmaps (presence and weight-1) over 240-number
words are scanned from the top with `ctz`, giving decreasing `n`. It falls
back to `std::sort` if:

- `p <= 5` (the events are then not coprime to 30), or
- the bitmap would need `>= N log2 N` words, or
- the word count exceeds `SIZE_MAX`.

---

## 5. Bounds, types and invariants (proofs to re-verify)

1. **PiTable range** (D's table has size `y`): every lookup in
   `DHyperbola.hpp` is `<= y`.
   - `pi[L]`, `pi[H]`, `pi[X/n]` (`< y`).
   - `pi[t/p]`, `pi[isqrt(t)]`.
   - `pi[min_r..max_r]` (`<= K/(min_n+1) < y/p`).
   - `pi[min_r2..max_r2]` (`< sqrt(y)`).
   - `pi[xpp/r]` (`<= max_n`), `pi[X/r^2]` (`< max_n`).
   - `pi[min_n/q]`, `pi[max_n/q]`, `pi[isqrt(max_n)]`.
   - `pi[y]` in `phi_minus_mertens(y, ...)`.
   - `pi[iroot<3>(xy)]` is `<= y` whenever `y >= iroot<3>(x)` (callers
     ensure `y > x^(1/3)`, except `x <= 3` with `y = 1`).
2. **64-bit results of `fast_div64`.** Every quotient is `< y` or
   `<= xy`:
   - `X/n` and `X/q` with `n, q > L`.
   - `K/max_n`, `X/max_n`.
   - `X/r^2` with `r > min_r2`.
   - `X/p^2 = x/p^3 <= x/y = xy` for primes that use the hyperbola method.
   - `X/y = L < y`.

   `fast_div(xp, high)`, `fast_div(xp, low1)` and `fast_div(xp, square)`
   return `T` and are clamped by `min(T, int64)`, which returns the
   `int64` type after comparing in `T` (see `include/min.hpp`).
3. **Products in `int64_t`:**
   - `r * r` with `r <= max_r2 < sqrt(y)`.
   - `q * primes[j] <= max_n`.
   - `square = p*p <= y`.
4. **Sieve-count order and range.** Every `count(t - low)` has
   `low <= t < high`:
   - For events, `n > min_n >= floor(X/high)` ⟹ `X/n < high`, and
     `n <= max_n <= floor(X/low)` ⟹ `X/n >= low`.
   - The order is `L <= H <= X/n` (decreasing `n`), followed by the
     kernel's leaves at positions `>= y`.
   - `test/gourdon/DHyperbola.cpp` asserts both properties.
5. **Signedness/modularity.** `T` is unsigned in all production paths
   (`UT`) and in the tests (`uint64_t`/`uint128_t`). Negative `int64`
   operands (`mertens_y`, `M(H)`, `-4b+2` terms) wrap consistently.
6. **Edge cases already tested:**
   - `x <= 3`, `y = 1` (`pi_gourdon`, D test);
   - `y = sqrt(x)`, incl. `x = 23^6 + 1` (`U = p^3` issue of PR #146,
     paper Section 6);
   - `floor(x/p^3) = y` (hard leaves only, `L = p^2`);
   - `L = H`, `L = p^2`;
   - prime fourth powers (`x = p^4`);
   - first composite `m` beyond `p^2`;
   - segment sizes 240 and 1680 in the dense test;
   - 32- and 64-bit prime vectors, 64- and 128-bit `T`.
7. **Known fragility removed by the restructuring.** An earlier helper
   computed `int64_t min_n = max(fast_div(xp, high), lower)`, which
   narrows a 128-bit value. It was unreachable for `x <= 10^31`
   (`LoadBalancerS2` segment sizes), but it is gone now. Look for similar
   patterns elsewhere.

---

## 6. Verification performed and how to reproduce

At `23b8549`, on Linux, with GCC 13:

- `-Wall -Wextra -pedantic -Werror` Release build with `BUILD_TESTS=ON`:
  - `ctest` 67/67;
  - `./primecount --test` passes;
  - `./primecount 1e18 -s` = `24739954287740860`.
- **Debug build** (`ENABLE_ASSERT`), as in CI:
  ```sh
  cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON \
        -DCMAKE_CXX_FLAGS="-march=x86-64 -Wall -Wextra -pedantic -Werror" \
        -DCMAKE_C_FLAGS="-Wall -Wextra -pedantic -Werror"
  cmake --build build-debug --parallel && (cd build-debug && ctest && ./primecount --test)
  ```
  All pass. The Debug `--test` takes about 3 minutes here.
- **Random sweep** (Appendix A.3): `primecount x --alpha=a -t{1,2,4}`
  versus `primecount x --deleglise-rivat`:
  - 400 cases with `10^5 <= x <= 3*10^14`;
  - 24 cases with `10^15 <= x <= 3*10^17`;
  - 400 cases with the Debug binary, `10 <= x <= 3*10^12`.

  `alpha` is random in `[1, x^(1/6)]`, so small `y` (`y < x^(2/5)`) is
  covered.
- **ARM:** compile-only checks with
  `clang++ --target=aarch64-linux-gnu -fsyntax-only -Werror` for
  `-march=armv8.2-a+sve` (native SVE), `-march=armv8-a` (NEON) and
  `-march=armv8-a -DENABLE_MULTIARCH_ARM_SVE`. They cover `D.cpp`,
  `AC.cpp`, `test/gourdon/{DBand,D_arm_sve,AC_arm_sve}.cpp`. Not executed
  on ARM hardware.
- **128-bit path consistency:** at `x = 10^20` (`T = int128_t`),
  `AC + D` equals `2919564469448030569` for both PR #146
  (AC `894394369769128121`, D `2025170099678902448`) and the extension
  (AC `866625209337419041`, D `2052939260110611528`).
- **Brute-force identity checks:**
  - Appendix A.1: classification and the cutoff `F`.
  - `test/gourdon/DHyperbola.cpp`: hyperbola identity against dense
    `Phi_p`/`M_p` tables, per segment, and the hyperbola method with easy
    leaves for every prime `(x/y)^(1/3) < p <= sqrt(y)` with
    `x <= 30000`.
  - `test/gourdon/D.cpp` `D_reference()`: whole-D brute force for small
    `x`.

After the C1 removal (`60db7d7d`), on the same machine:

- Release and Debug builds as above: `ctest` 67/67, `--test` passes and
  `1e18` is correct.
- `--AC` and `--D` give the same values as at `23b8549` for 300 random
  `(x, alpha, threads)` cases with `x <= 3*10^14`, for `x = 2^64 + 1`
  (default `alpha`) and for `x = 3*10^19 + 7` (`alpha = 20`).
- Random sweep (Appendix A.3): 200 cases with `10^5 <= x <= 3*10^14`.
- The ARM compile-only checks above, and `AC.cpp` on x86 with and
  without `ENABLE_LIBDIVIDE`.

After the AC primes reduction (`bfb74102`):

- Release and Debug builds: `ctest` 67/67 (including the new
  `prev_prime()` test in `test/PiTable.cpp`), `--test` passes and `1e18`
  is correct.
- `--AC` gives the same values as at `60db7d7d` for 400 random
  `(x, alpha, threads)` cases with `x <= 3*10^16`, 400 cases with the
  Debug binary (`Vector` index asserts) and `x <= 10^12`, and for
  `x = 2^64 + 1`, `3*10^19 + 7` (`alpha = 20`) and `10^20 + 3`
  (`alpha = 3`).
- Random sweep (Appendix A.3): 200 cases with `10^5 <= x <= 3*10^14`.
- ARM compile-only checks as above, and x86 builds of `AC.cpp` with and
  without `ENABLE_LIBDIVIDE` and with `DISABLE_INT128`.

After the switch to the hyperbola terminology (`22f2be25`, no
functional change): Release and Debug builds pass `ctest` 67/67 (the
test is now called `DHyperbola`) and `--test`, `1e18` is correct, and
the ARM compile-only checks pass.

Not yet done (recommended):
- run on real ARM (SVE and NEON);
- MSVC and macOS builds;
- a larger-`x` 128-bit run against a known `pi(x)`, e.g.
  `primecount 1e21` = `21127269486018731928`.

---

## 7. Performance data

Census of D's sieve queries, from a separate counting tool (not in the
repo). It models the extension as of `b9ef7c5`, with the primes
`floor(x/p^3) = y` counted as having easy leaves, which differs
negligibly from the current code.

| x | y (default) | master | PR #146 | extension | removed (PR) | removed (ext.) | C1 easy leaves absorbed |
|---|---|---:|---:|---:|---:|---:|---:|
| 1e16 | 4,496,510 | 111,012,964 | 92,108,987 | 88,463,591 | 17.03% | 20.31% | 16,540,992 |
| 1e18 | 32,222,000 | 1,524,253,981 | 1,178,822,430 | 1,116,680,620 | 22.66% | 26.74% | 332,530,670 |
| 1e20 | 245,535,363 | 22,178,370,075 | 15,209,216,839 | 14,068,992,941 | 31.42% | 36.56% | 7,693,203,942 |
| 1e22 | 1,885,604,250 | 329,139,642,154 | 196,456,547,550 | 176,703,402,636 | 40.31% | 46.31% | 168,430,469,271 |
| 1e24 | 14,071,600,000 | 4,851,837,541,432 | 2,545,923,003,175 | 2,228,535,552,416 | 47.53% | 54.07% | 3,323,611,530,258 |

At `1e20` the `E_p` prime terms (PiTable lookups) grow from 2.37e9
(PR) to 4.73e9, while C1's 7.69e9 PiTable lookups disappear from AC.

Timings come from a 4-core Xeon VM with 4 threads. This machine is noisy
(±5–10%), so re-measure on real hardware.

- At `x = 10^20`:
  - D: 83.2 s (PR #146) vs 83.5 s (extension);
  - AC: about 65 s vs 57–60 s.
- At `x = 10^19`, the internal counts are:
  - `hyperbola_sum()` calls: 3.2M (pairs of prime and segment), of which
    111k with hard leaves only and 16k with easy leaves;
  - events: 78M;
  - prime-`n` terms in segments: 818M;
  - remaining prime leaves of the primes with easy leaves: 2.1M;
  - the easy-terms loop is 1.07 CPU-s out of about 135 CPU-s for D (0.8%).
- The structural refactoring (`f257c8c`) is performance-neutral. D at
  1e19 took 36.4 s vs 36.9 s (mean of 4). Callgrind instructions at
  1e16 changed by −0.8% (`-O2`) and +0.5% (`-O3`), i.e. compiler
  noise.

---

## 8. Design decisions and rejected alternatives

- **The leaves `(X/y, y]` instead of the natural region `(p^2, x/p^3]`.**
  - It reaches the frontier: for every prime that uses the hyperbola
    method, all its leaves with position `< y` are in `(X/y, y]`, so all
    their prime-leaf sieve queries cancel.
  - It covers the primes `x^(1/5) < p <= sqrt(y)`, which did not use the
    hyperbola method in PR #146 (`H >= x/p^3` there).
  - It turns C1's per-leaf `pi` lookups into folded `E_p` terms.
- **Primes with `floor(x/p^3) = y` have hard leaves only.** They have no
  easy leaves, and `L = p^2`.
- **Easy hyperbola terms are computed once per prime**, not inside the
  segment containing `p^2`.
  - Threads with limit `< p^2` skip the prime.
  - Putting the terms into the segment containing `p^2` breaks edge
    cases: `p^4 <= x < p^3 (p+1)` (`goto next_segment`) and `y = p^2`
    with `low = p^2` (cut-off).
  - It would also put millions of `pi` lookups on a single thread. The
    `RelaxedAtomic` loop instead balances them across threads.
- **The remaining prime leaves of primes with easy leaves go through the
  factor table** rather than a dedicated prime loop. This gives one code
  path, no `bool` return and no duplicated `phi[b]` / `cross_off_count`.
  It is performance-neutral (2.1M queries at 1e19).
- **Rejected:**
  - A `DHyperbola` class (fewer parameters, but boilerplate cancels the
    savings).
  - libdivide for the easy-terms loop (only 0.8% of D).
  - Caching per-prime values across segments (negligible).
  - Replacing the sorter by a factor-table scan: about 10x more entries
    scanned than there are events.
- **AC-only alternative** (branch `claude/busy-cori-gwffdl`, commit
  `a819adb`): keeps D unchanged and applies Gourdon's inversion to C1's
  clustered leaves. It is superseded by the full extension, and moot now
  that C1 is removed.

---

## 9. Review checklist for Codex

Verify each item and report issues with a concrete counterexample.

1. Section 3.1:
   - `E_p` for `p^2 <= t < p^3`;
   - `Phi_p`, `M_p`, `E_p` for `p <= t < p^2`;
   - check numerically against dense tables.
2. Section 3.2: the hyperbola identity with `L < p^2` (easy leaves),
   including the weight-1 event `n = p^2`.
3. Section 3.3:
   - every `m` in `(L, y]` is a leaf;
   - for primes with easy leaves, the leaves `L < m <= y` contain exactly
     C1's leaves of `p`;
   - D's primes with easy leaves are exactly C1's primes
     `(x/y)^(1/3) < p <= sqrt(y)`.
4. Section 3.4: the three integer equivalences and the monotonicity. Also
   check `F <= iroot<3>(xy)` (Section 10.1), on which the easy-terms loop
   of `D.cpp` relies.
5. Section 3.5: each hyperbola term is computed exactly once, either in
   `easy_hyperbola_terms()` or in exactly one segment by
   `hyperbola_sum()`, across all `LoadBalancerS2` work units (single- and
   multi-threaded segment sizing), `goto next_segment` and the cut-off.
6. Section 3.7: Gourdon's inversion with floors in
   `sum_phi_minus_mertens()`, including `min_r2`/`max_r2` via `isqrt`.
7. Section 5: PiTable ranges, 64-bit quotient bounds, `int64_t` products,
   narrowing in `min`/`max`, unsigned arithmetic, implicit signed conversions.
8. Kernels:
   - query order (`L`, `H`, events, then factor-table leaves);
   - `max_m` lowering;
   - the cut-off before `hyperbola_sum()` and its `low >= y` condition;
   - the 2-prime loop `max3(xp_high, y / prime, prime)`;
   - all five kernels are identical in these parts.
9. `D.cpp`:
   - `mertens_y` validity for every prime that reads it;
   - the easy-terms loop bounds equal the primes with easy leaves;
   - `PiTable`/primes sized to `y`;
   - `RelaxedAtomic` usage.
10. `SemiprimeEventSorter`: correctness of the bitmap path (mod-240
    wheel, `p > 5`), the fallback condition, and 32-bit `size_t`
    behaviour.
11. Tests: whether `test/gourdon/DHyperbola.cpp` and `D.cpp` cover primes
    with easy leaves and `H >= p^2` or `H < p^2`, primes with hard leaves
    only, 128-bit `T` with
    `xp > 2^64` (currently only through `--test`/`D_128bit` with small
    values), and `y < x^(2/5)`. Add tests where coverage is missing.
12. Conventions: `AGENTS.md` (C++14, internal `min`/`max`, include order,
    comments, copyright year, ChangeLog). Run the Release and Debug
    verification (Section 6).

---

## 10. C1 removal: C1 was always empty

### 10.1 The claim

**Correction.** Earlier versions of this section claimed that C1 is
empty only if `y >= x^(2/5)`, and proposed a lower bound for `y`. That
was wrong. The argument compared `band_min_prime()` only with C1's lower
bound `(x/y)^(1/3)` on `p` and ignored its upper bound `p <= sqrt(y)`, and
the check script only sampled `y >= x^(2/5)`. Appendix A.2 now samples
every valid `y`. `C1-Removal-Analysis.md` found the error; this section
follows its proof.

**Claim:** for all integers `1 <= y <= N`, the cutoff `F` of Section 3.4
(formerly `band_min_prime(N, y)`) satisfies:

```
F = min(max(iroot<3>(y), N/y), isqrt(y)) <= iroot<3>(N)
```

Every valid Gourdon setup has `y^2 <= x`, hence `y <= N = floor(x/y) = xy`.
This includes tiny `x` (`x = 2, 3` with `y = 1`, `x = 8` with `y = 2`).

**Proof** (`C1-Removal-Analysis.md`, Section 4.1). Let `R = iroot<3>(N)`.

- `iroot<3>(y) <= R`, because `y <= N`.
- If `isqrt(y) <= R`, the outer `min()` gives the claim.
- Otherwise `y >= (R+1)^2`. Since `N < (R+1)^3`, we get
  `N/y < (R+1)^3/(R+1)^2 = R+1`, so `floor(N/y) <= R`.

**Consequences.**

- AC: `max_c1 <= pi_root3_xy < min_c1`, so the C1 loop was never
  entered (Section 4.5).
- Every C1 prime uses the hyperbola method and has easy leaves. A C1
  prime has `p^2 <= y` and `x < p^3*y` (⟺ `p > iroot<3>(xy)`). Together
  with `y^2 <= x`:
  - `x < p^3*y <= p*y^2`, so `H < y`;
  - `y^2 <= x < p^3*y`, so `y < p^3`;
  - `x < p^3*y`, so `L = floor(x/(p*y)) < p^2`;
  - `p^4 <= y^2 <= x`, so `p` is below D's cap `x^(1/4)`.

  The leaves `L < m <= y` include all of C1's leaves of `p` (Section 3.3).
- `y > x^(2/5)` is necessary for C1 primes and for primes that use the
  hyperbola method (Section 3.4), but it is not a condition for removing
  C1. For `y <= x^(2/5)`, C1's range is empty, already in master. For
  `y > x^(2/5)`, D computes the easy leaves of every C1 prime with the
  hyperbola method.

**Verification.**

- `C1-Removal-Analysis.md`, Section 10: exhaustive and random checks of
  the claim in C++ (against the then-production `band_min_prime()`),
  JavaScript and C#.
- Python with exact integers:
  - the claim for all `1 <= y <= N <= 2500` and for 200,000 random
    `x < 10^31` (3,326,250 cases);
  - 200,000 random primes, including 24,736 C1 primes. Each C1 prime
    uses the hyperbola method with easy leaves, and has `p <= x^(1/4)`
    and `y^5 > x^2`.
- A build with `abort()` in the C1 loop body of all three AC kernels
  never aborted:
  - `ctest` 67/67, including the AC tests with small `y`;
  - `primecount --test`;
  - 400 random `(x, alpha, threads)` cases with `x <= 3*10^14` against
    `--deleglise-rivat`;
  - `--AC --alpha={1,2,3.5}` at `x = 10^6, 10^9, 10^12, 10^15`.

### 10.2 Current defaults

`y = floor(floor(x^(1/3)) * alpha)`, with `alpha = get_alpha_gourdon(x)`,
clamped to `[x^(1/3) + 1, sqrt(x)]`. Primes `p > k` that use the
hyperbola method (`k = 8` for all `x` below):

| x | default y | y / x^(2/5) | hard leaves only | also easy leaves (= C1 primes) |
|---|---:|---:|---:|---:|
| 1e6 | 207 | 0.824 | 0 | 0 |
| 1e9 | 2,619 | 0.658 | 0 | 0 |
| 1e11 | 13,825 | 0.550 | 0 | 0 |
| 4e11 | 43,824 | 1.002 | 0 | 0 |
| 4.1e11 | 44,545 | 1.009 | 0 | 1 |
| 1e12 | 76,500 | 1.212 | 12 | 7 |
| 1e20 | 245,535,363 | 2.455 | 680 | 888 |
| 1e31 | 11,631,189,649,596 | 4.630 | 67,681 | 169,271 |

Below about `x = 4*10^11`, D does not use the hyperbola method at default
parameters, and C1 had no primes in master either. Above it, the C1
primes are D's primes with easy leaves. No change of `alpha` or `y` is
needed.

### 10.3 Implementation

`1aeaa05f` and `60db7d7d` implement `C1-Removal-Analysis.md`, Sections
9.1 and 9.2. `alpha`, `y`, `get_alpha_gourdon()` and the seven places that
compute Gourdon's `y` (`pi_gourdon.cpp` and `main.cpp`) are unchanged
(`C1-Removal-Analysis.md`, Section 8).

1. **C1 deleted** (`1aeaa05f`):
   - `C1()` in `AC_default.hpp`, `C1_libdivide()` in `AC_libdivide.hpp`,
     `C1_arm_sve()` in `AC_arm_sve.hpp`;
   - the `if (low < y) { ... C1 loop ... }` blocks, `max_c1` and `sqrty`
     in the three `AC_OpenMP_*` functions;
   - `#include "DBand.hpp"` in AC.

   `LoadBalancerAC`, `SegmentedPiTable` and the segments below `y` are
   kept: A and C2 still have work there. The thread-reuse comments list
   only C2 and A, and the C2 doc comments say that D computes the first
   part of C.
2. **Tests.** All AC and D test values are unchanged, including those
   with small `y`, e.g. `{1000000, 101}`, `{1000000, 207}`,
   `{10000000, 485}`. Their C1 range was empty, already in master.
3. **D** (`60db7d7d`): the easy-terms loop starts at
   `max(k, pi_root3_xy) + 1` instead of `max3(k, pi_root3_xy, pi_band) + 1`,
   and `band_min_prime()` is removed. `test/gourdon/DBand.cpp` (now
   `DHyperbola.cpp`) checks instead that every prime
   `(x/y)^(1/3) < p <= sqrt(y)` uses the hyperbola method with easy
   leaves, exhaustively for `x <= 30000`.
4. The kernel cut-off's `low >= y` condition is unchanged and still
   required (Section 4.3).
5. `ChangeLog`: "AC_*.hpp: Remove the C1 formula, D computes its easy
   leaves using the hyperbola method." The tex paper is not updated yet
   (Section 11).
6. **Verification:** Section 6.

---

## 11. Open items and refinement ideas

- **Cheaper test in `hyperbola_sum()`.** About 96% of calls are for
  primes that don't use the hyperbola method. Test `xpy >= y`
  (⟺ `isqrt(xp) >= y`, Section 3.4) first, using
  `xpy = min(fast_div(xp, y), y)` to stay 64-bit safe, and compute
  `isqrt(xp)` only for the primes that use it.
- **`easy_hyperbola_terms()`** recomputes `L`, `H` and `X/p^2` per prime;
  this is harmless. Its prime loop could use libdivide; probably not
  worth it (0.8% of D).
- **Test coverage.** Real `x > 2^64` with 128-bit `T` for primes that use
  the hyperbola method (the dense test only uses small `xp` with
  `uint128_t`).
- **Remove the `p` parameter** of `SemiprimeEventSorter::sort` if
  `p > 5` is always guaranteed (k >= 3?). Check `PhiTiny::get_k()` for
  small `x`; currently `p <= 5` occurs for tiny `x` and needs the
  `std::sort` fallback.
- **Revisit the paper** (§6–8) to use the terms of Section 2 (the
  hyperbola method, primes with hard leaves only or also easy leaves)
  instead of the natural region, the y-region, band and row, and to
  describe the C1 removal.

---

## Appendix A. Verification scripts

### A.1 Hyperbola test versus the cutoff `F` (Python)

```python
import math, random

def iroot3(n):
    r = round(n ** (1 / 3))
    while r ** 3 > n: r -= 1
    while (r + 1) ** 3 <= n: r += 1
    return r

def check(x, y, p):
    xp = x // p
    sqrt_xp = math.isqrt(xp)
    uses_hyperbola = sqrt_xp < y and p * p > y // p    # hyperbola_sum()
    xy = x // y
    bound = min(max(iroot3(y), xy // y), math.isqrt(y))  # F, Section 3.4
    assert uses_hyperbola == (p > bound), (x, y, p)
    if uses_hyperbola:
        assert ((xp // y) < p * p) == (p > iroot3(xy)), (x, y, p)  # easy leaves

rng = random.Random(1)
for _ in range(200000):
    x = max(int(10 ** rng.uniform(3, 31)) + rng.randrange(-5, 6), 8)
    lo, hi = iroot3(x) + 1, math.isqrt(x)
    if lo > hi: continue
    y = rng.randint(lo, hi)
    check(x, y, rng.randint(2, math.isqrt(y)))
for p in range(2, 3000):
    for x in {p**4 - 1, p**4, p**4 + 1, p**5 - 1, p**5, p**5 + 1, p**6 - 1, p**6, p**6 + 1}:
        for y in {p*p - 1, p*p, p*p + 1, p**3 - 1, p**3, p**3 + 1, math.isqrt(x), iroot3(x) + 1}:
            if iroot3(x) < y <= math.isqrt(x) and p * p <= y:
                check(x, y, p)
print("OK")
```

### A.2 `F <= iroot<3>(xy)` for every valid `y` (Python)

Checks the claim of Section 10.1 for random `y` in the whole valid range,
below and above `x^(2/5)`.

```python
import math, random
from decimal import Decimal, getcontext
getcontext().prec = 80

def iroot(n, k):
    r = int(Decimal(n) ** (Decimal(1) / k))
    while r ** k > n: r -= 1
    while (r + 1) ** k <= n: r += 1
    return r

rng = random.Random(3)
for _ in range(300000):
    x = max(int(10 ** rng.uniform(2, 31)) + rng.randrange(-3, 4), 2)
    ymin = min(iroot(x, 3) + 1, math.isqrt(x))
    y = rng.randint(ymin, math.isqrt(x))
    xy = x // y
    F = min(max(iroot(y, 3), xy // y), math.isqrt(y))
    assert F <= iroot(xy, 3), (x, y)
print("OK")
```

### A.3 Random sweep against Deléglise–Rivat (Python)

```python
import random, subprocess, sys
NEW = sys.argv[1]            # primecount binary under test
REF = sys.argv[2]            # reference binary (any build), uses --deleglise-rivat
def run(b, args): return subprocess.run([b] + args, capture_output=True, text=True, check=True).stdout.split()[0]
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
for _ in range(400):
    x = int(10 ** rng.uniform(5, 14.5)) + rng.randrange(1000)
    alpha = round(rng.uniform(1, max(1.0, x ** (1 / 6))), 3)
    got = run(NEW, [str(x), f"--alpha={alpha}", f"-t{rng.choice([1, 2, 4])}"])
    exp = run(REF, [str(x), "--deleglise-rivat", "-t4"])
    assert got == exp, (x, alpha, got, exp)
print("OK")
```
