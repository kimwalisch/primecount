# Removing C1 without changing Gourdon's tuning parameters

Analysis for independent verification by Claude, 2026-10-22

Source revision: `8aaad815efb1f88b29753354fffe9495cba1c4ff`, branch
`claude/d-extended-band`. Function names and formulas below identify the relevant
code; line numbers refer to this revision and may move later.

## 1. Conclusion and scope

**The C1 implementation in AC is already unreachable for every currently valid
choice of `y`. Removing it does not require imposing `y >= x^(2/5)`.**

The full mixed-band extension on this branch already transfers every possible
C1 row into D. For smaller `y`, the original C1 row interval is empty even before
the transfer. This is stronger than the sufficient condition proved in
[Hard-Special-Leaves-Band.md, Section 10](Hard-Special-Leaves-Band.md).

The recommendation is to delete the unreachable C1 implementations and their
dispatch blocks, preserve the existing tuning and test expectations, and correct
the design document. No new alpha getter, `y` clamp, parameter object, or rational
power computation is needed for this removal.

The user's follow-up question was whether retaining `y < x^(2/5)` requires D to
allocate primes beyond `y` or a PiTable covering arguments greater than `y`.
The answer is also no, subject to the existing algorithm's lower bound near
`x^(1/3)`. Sections 6 and 7 explain that separately from the C1 proof.

This document is an analysis and implementation proposal. It does not implement
the removal, re-prove the entire band identity, or replace a complete review of
the new D algorithm. In particular, removing unreachable C1 code preserves the
current program's behavior; the correctness of the existing D transformation is
a separate obligation, supported by its derivation and existing reference tests.

Claude should independently challenge the assumptions, proofs, and code mapping
below before implementation starts. The user has requested a document for that
review and has not yet authorized the algorithm changes.

## 2. Sources and current code structure

Relevant mathematical background:

- [Hard-Special-Leaves-Band.md](Hard-Special-Leaves-Band.md), especially
  Sections 3.3, 3.4, 4.3, 4.5, 5, and 10. It describes this branch's full
  mixed-band extension.
- [Hard_Special_Leaves_Cancellation_2.tex](Hard_Special_Leaves_Cancellation_2.tex),
  especially the band rearrangement and the distinction between the natural
  region and the y-region. It describes the earlier algorithm; its region
  definitions must not be substituted for the full extension now in the code.
- [Easy-Special-Leaves.pdf](Easy-Special-Leaves.pdf), for AC, segmentation,
  and the reflected C2 contribution.
- [Hard-Special-Leaves.pdf](Hard-Special-Leaves.pdf) and
  [Hard-Special-Leaves-SIMD-Filtering.pdf](Hard-Special-Leaves-SIMD-Filtering.pdf),
  for the underlying segmented hard-leaf calculation.
- The original Gourdon paper is available in the sibling `../math_papers`
  directory relative to the repository root. Its general `y,z` formulation
  should be distinguished from this branch, which uses `z = y`.

Principal implementation references:

| Location | Role |
| --- | --- |
| [AC_default.hpp](../src/gourdon/AC_default.hpp), `C1()` and `AC_OpenMP_default()` | Scalar C1 implementation and row bounds; relevant bounds start near lines 310 and 364 |
| [AC_libdivide.hpp](../src/gourdon/AC_libdivide.hpp), `C1_libdivide()` and `AC_OpenMP_libdivide()` | Same row bounds, with scalar fallback for wide divisions |
| [AC_arm_sve.hpp](../src/gourdon/AC_arm_sve.hpp), `C1_arm_sve()` and `AC_OpenMP_arm_sve()` | Same row bounds for ARM SVE |
| [DBand.hpp](../src/gourdon/DBand.hpp), `band_min_prime()` | Integer band cutoff, near line 40 |
| [DBand.hpp](../src/gourdon/DBand.hpp), `easy_band_terms()` and `band_row()` | Split between easy band terms and terms requiring sieve queries |
| [D.cpp](../src/gourdon/D.cpp), `D_OpenMP()` | Prefix preparation, row scheduling, and easy-band loop, near lines 175-216 |
| [D_default.hpp](../src/gourdon/D_default.hpp), `D_thread_default()` | Representative D kernel; the corresponding SIMD kernels have the same mathematical bounds |
| [pi_gourdon.cpp](../src/gourdon/pi_gourdon.cpp) | 64-bit and 128-bit full-algorithm parameter construction |
| [main.cpp](../src/app/main.cpp) | Five CLI partial-formula parameter construction paths |
| [util.cpp](../src/util.cpp), `get_alpha_gourdon()` and `get_alpha()` | Requested/default tuning selection and effective-alpha reporting |
| [phi_vector.cpp](../src/phi_vector.cpp) | Initialization of partial sieve counts, including bounded PiTable shortcuts |

## 3. Parameters and exact integer notation

The current implementation has one Gourdon split parameter: `z = y`. For
ordinary inputs the existing construction enforces

```text
iroot<3>(x) + 1 <= y <= isqrt(x).
```

Thus `x^(1/3) < y <= sqrt(x)`. The upper clamp takes precedence for tiny inputs:
for example, `x = 2,3` use `y = 1`, and `x = 8` uses `y = 2`.
CLI partial formulas also accept `x = 1`; full prime counting returns early
for `x < 2`.

All divisions below are floor divisions of positive integers. Define:

```text
N = floor(x/y)                  // named xy in the code
R = floor(cuberoot(N))          // iroot<3>(xy)
S = floor(sqrt(y))              // isqrt(y)
Q = floor(cuberoot(y))          // iroot<3>(y)
A = floor(N/y)                  // floor(x/y^2)
F = min(max(Q, A), S)           // band_min_prime(N, y)
```

For a row prime `p = primes[b]`, define:

```text
X = floor(x/p)
H = floor(sqrt(X))
L = floor(X/y)
```

The integer identity `floor(floor(x/p)/y) = floor(x/(p*y))` is used throughout.
Statements such as `x < p^3*y` below are mathematical inequalities, not a
recommendation to evaluate their products in fixed-width production arithmetic.

There are two different assumption sets:

1. **C1 dead-code theorem:** only `x >= 1`, `y >= 1`, and `y^2 <= x` are needed
   for the cutoff inequality proved in Section 4.
2. **Validity of the complete algorithm and its allocations:** retain all
   existing input restrictions, the lower `y` bound above, and the 128-bit
   input guard `x <= (int128_t(y) << 62)` where currently enforced.

The first theorem does not establish support for arbitrary `y` below the
existing cube-root bound. Nor does it allow removing the large-input guard.
The supported maximum input remains `10^31`.

## 4. Exact proof that the current C1 loops are empty

### 4.1 The cutoff inequality

Claim:

```text
F = min(max(iroot<3>(y), N/y), isqrt(y)) <= iroot<3>(N) = R
```

for every pair of positive integers `y <= N`.

The current upper bound `y^2 <= x` gives `y <= floor(x/y) = N`, so the claim
applies to every currently valid Gourdon setup, including the tiny cases.

Proof:

1. Since `y <= N`, monotonicity of the integer cube root gives `Q <= R`.
2. If `S <= R`, the outer `min(..., S)` immediately gives `F <= R`.
3. Otherwise `S >= R+1`, so `y >= (R+1)^2`. By the definition of `R`,
   `N < (R+1)^3`. Therefore

   ```text
   N/y < (R+1)^3 / (R+1)^2 = R+1,
   ```

   and the integer quotient `A = N/y` satisfies `A <= R`.
   Both `Q` and `A` are then at most `R`, giving `F <= R` again.

No floating-point approximation or assumption about `x^(2/5)` is involved.

### 4.2 Mapping the inequality to AC

Each of the three AC backends computes:

```cpp
int64_t pi_root3_xy = pi[iroot<3>(xy)];
int64_t max_c1 = pi[band_min_prime(xy, y)];
```

The per-segment C1 starting index includes:

```cpp
int64_t min_c1 = max(k, pi_root3_xy);
// Additional segment bounds can only increase min_c1
min_c1 = max3(min_c1, pi_sqrt_low, pi[min_c1_prime]) + 1;
```

Monotonicity of prime counts and Section 4.1 give:

```text
max_c1 = pi(F) <= pi(R) < min_c1.
```

Thus the `for (b = min_c1; b <= max_c1; b++)` body is never entered.
The `low < y` condition and other segment restrictions only remove more work.
They cannot invalidate this argument.

This is a direct proof about the existing source. Deleting these helpers and
loops requires no change to the values of AC or D and no regeneration of
expected results.

## 5. Why D covers every original C1 row

The original C1 candidate range, before this branch's band exclusion, is:

```text
max(k, pi(R)) < b <= pi(S).
```

For any candidate prime `p`, this gives:

```text
p^2 <= y
p^3 > floor(x/y), equivalently x < p^3*y.
```

The equivalence is exact: an integer `p^3` exceeding `floor(x/y)` also exceeds
the real quotient `x/y`.

Combining these conditions with `y^2 <= x` gives:

```text
x < p^3*y <= p*y^2       => floor(sqrt(floor(x/p))) < y
y^2 <= x < p^3*y        => y < p^3.
```

These are precisely the two conditions checked by `band_row()`, in the range
`p^2 <= y`. Moreover `x < p^3*y` implies `L = floor(x/(p*y)) < p^2`, so this
is a mixed band, which includes the row's easy leaves.

D's row cap also loses nothing: `p^4 <= y^2 <= x`, giving
`p <= iroot<4>(x)`. The common `k` bound excludes the same initial rows.
The equality case at a fourth power is included by that cap; it introduces
no missing C1 row.

For completeness, C1's easy leaves have

```text
max(x/p^3, y/p) < m <= y, with lpf(m) > p and mu(m) != 0.
```

Since `p^2 <= y` and `y^2 <= x`, we have `p^2*y <= x`, hence
`x/p^3 >= y/p`. Also `p <= y`, so `L <= floor(x/p^3)`.
Every such easy leaf therefore lies in `(L, y]` and is included by the mixed
band. At its query position `t = floor(x/(p*m)) < p^2`, its value is
`-mu(m) * (pi(t) - b + 2)`, the same contribution as in C1.

This describes the transfer already implemented on the branch. It does not
propose doing the transfer a second time.

### 5.1 What the `x^(2/5)` threshold actually means

An original C1 candidate would require

```text
x < p^3*y <= y^(5/2).
```

Therefore original C1 candidates can occur only when `y > x^(2/5)`.
Below or at that threshold, the original C1 interval is already empty.
Above it, every possible original C1 row belongs to D's mixed band by the
argument above. Primality and `k` can make the interval empty even above
the threshold.

The same threshold is necessary for *any* cancellation band on this branch,
including a hard band. An active band requires `H < y` and `p^2 <= y`.
The first condition is equivalent to `x < p*y^2`, so

```text
x < p*y^2 <= y^(5/2).
```

Consequently, for `y <= x^(2/5)`, all calls to `band_row()` return at the
initial guard, and the `easy_band_terms()` row interval is also empty.
D computes the leaves with its ordinary segmented-sieve path.

### 5.2 Concrete contradictions to the design document

The following are explicit low-`y` examples named in Section 10:

| `x` | `y` | `R = iroot<3>(x/y)` | `S = isqrt(y)` | Original C1 range before `k` |
| --- | --- | --- | --- | --- |
| 1,000,000 | 101 | 21 | 10 | `21 < p <= 10`, empty |
| 1,000,000 | 207 | 16 | 14 | `16 < p <= 14`, empty |
| 10,000,000 | 485 | 27 | 22 | `27 < p <= 22`, empty |
| 1,000,000 | 1,000 | 10 | 31 | Potential rows, all transferred into mixed bands |

Thus Section 10.3's assertion that deleting C1 makes the first three AC cases
incomplete is incorrect. Their expected values should be preserved.

## 6. Why small `y` does not require a larger PiTable or prime array

The D entry points currently construct:

```cpp
PiTable pi(y, threads);
auto primes = pi.get_primes<uint32_t>(y, threads);
```

The 128-bit path makes the same choice of limit, using a wider prime element
type when necessary. The prime vector has `pi(y)+1` entries because index zero
is reserved. It does not contain an extra successor prime beyond `y`.

### 6.1 Setup and ordinary D kernels

For ordinary supported inputs, `y > x^(1/3)` ensures:

```text
iroot<4>(x) <= y
iroot<3>(floor(x/y)) <= y
isqrt(y) <= y.
```

For the second inequality, `x/y < y^2`, so the cube root of `x/y` is less
than `y^(2/3) <= y`. These bounds cover D's setup PiTable accesses and row
prime indices. The tiny exceptions also satisfy the displayed integer bounds
directly, even where the real cube-root lower bound fails.

In each D kernel:

- Row primes are capped by `iroot<4>(x)`.
- Factor-table leaves are clamped to `m <= y`.
- The ordinary two-prime loop sets `l = pi[max_m]`, with `max_m <= y`.
  Thus its `primes[l]` accesses stay within the allocated vector. An empty
  range can use `primes[0] = 0`; the code does not need a successor sentinel.

The query position `t = floor(x/(p*m))` can exceed `y` when `y` is small.
That is expected. D evaluates its partial-sieve contribution as

```cpp
int64_t count = sieve.count(t - low);
int64_t phi_t = phi[b] + count;
```

using the sieve segment containing `t`. It does not require `pi[t]`.
The segmented-sieve distance extends to `x/y`, independently of the PiTable's
limit `y`. These are different data structures with different ranges.

For `y <= x^(2/5)`, Section 5.1 establishes that the new band transformation
is inactive. Therefore there is no band-specific attempt to replace an
above-`y` sieve query by a PiTable lookup.

### 6.2 `phi_vector()` with a starting position above `y`

The D kernels initialize `phi_vector(low, max_b, primes, pi)` even when
`low > y`. This does not mean that they index `pi[low]` unconditionally.

`PhiCache::is_pix()` in `phi_vector.cpp` first checks
`x < pi_.size()`. The π shortcut is used only inside the supplied table's
range; otherwise the recursive/cached partial-sieve calculation is used.

The outer function also has a direct `pi[x]` lookup in the branch
`if (primes[a] > x) a = pi[x];`. This branch implies
`x < primes[a] <= y`, so that lookup is within the supplied table too.

The `primes_[a+1]` lookahead inside this predicate also fits. In
`phi_vector()`, the initial requested row count is at most `pi(iroot<4>(x))`.
Its recursive calls pass `i-2`, with `i` at most that row count. Thus the
lookahead is within the original supplied vector. The direct `primes[a]`
access in the outer function is within that row count as well.

The loop that enters its π-only continuation does so after finding a value
for which `is_pix()` succeeds; subsequent quotients decrease while the
relevant prime threshold increases. This is part of the existing
`phi_vector()` algorithm, not a new requirement introduced by C1 removal.

## 7. Bounds on active-band lookups, including reflection

Although no bands run below the threshold, active-band table safety can also
be established locally without imposing a new global `y` clamp.

### 7.1 Direct band arguments and semiprime events

An active row satisfies `H < y`. With `L = floor(X/y)`, we have `L <= H < y`.
The summation bounds satisfy:

```text
L <= min_n < max_n <= H < y.
```

Every actual summation value or event has `n > L`, so

```text
floor(X/n) < y.
```

Thus the direct `pi[L]`, `pi[H]`, and `pi[X/n]` lookups fit. The function
`phi_minus_mertens(t, p, b, pi)` accesses only `t`, `floor(t/p)`, and
`isqrt(t)`. Its band calls have `t < y`, and the prefix preparation calls
it with `t = y`.

Semiprime generation uses prime indices bounded by
`pi[isqrt(max_n)]` and `pi[max_n/q]`. Both arguments are at most `y`.
Every enumerated product `q*r` is at most `max_n < y`; neither prime factor
can exceed `y`.

### 7.2 Reflected prime ranges and quotient lookups

In `sum_phi_minus_mertens()`, let `K = floor(X/p)` and
`u = min_n`, `v = max_n`. The reflected endpoints are:

```text
min_r  = floor(K/v)
max_r  = floor(K/(u+1))
min_r2 = floor(sqrt(floor(X/v)))
max_r2 = floor(sqrt(floor(X/(u+1))))
```

Since `u+1 > L`, `floor(X/(u+1)) < y`. Therefore `max_r < y/p` and
`max_r2 < sqrt(y)`. The endpoint PiTable lookups and all prime indices
inside these reflected ranges fit in the existing allocations.

The inverse quotient accesses also fit, which is a separate check:

- A reflected prime `r > floor(K/v)` satisfies `floor(K/r) < v < y`.
- A reflected prime `r > floor(sqrt(floor(X/v)))` satisfies
  `floor(X/r^2) < v < y`.

The second statement follows because `r^2 > floor(X/v)` implies
`r^2 > X/v`. No assumption about the distribution of primes is needed.

The ordinary, non-reflected sum uses `n > L`, already bounded in Section 7.1.
Its four-way unrolling checks that all four indices are at most the last
index; it requires no prime-array padding beyond `pi(y)`.

### 7.3 Easy band terms and prefix validity

`easy_band_terms()` is called only for mixed-band rows. The preceding proofs
give `L < p^2 <= y` and `H < y`, so its direct endpoint lookups fit.
Its summation arguments exceed `floor(X/p^2)` and `L`; every resulting
query position is below `y`.

One distinction is important for an independent review: D initializes
`mertens_y[b]` for a larger prefix of row indices than it actually consumes
in the band transformation. The closed formula for
`phi_minus_mertens(y, p, b, pi)` is mathematically valid only for
`p^2 <= y < p^3`. Entries computed outside that domain must not be described
as correct restricted Mertens values. The current callers consume them only
in that domain:

- `band_row()` is called from the `p^2 <= y` part of the kernel and returns
  before using its prefix argument unless `y < p^3` and `H < y`.
- `easy_band_terms()` requires `p > R`, which gives `p^3 > N >= y`, as
  well as `p^2 <= y` from its upper row bound.

Computing the unused entries is table-safe because their lookup argument is
still `y`; their mathematical validity and their memory bounds are distinct
questions. Restricting this precomputation would be a separate optimization.

## 8. Alpha API analysis

Seven current setup paths obtain `alpha`, convert it to `y`, then clamp `y`:
two in `pi_gourdon.cpp` and five partial-formula CLI paths in `main.cpp`.

`get_alpha_gourdon()` currently selects the default or requested tuning,
applies the `--double-check` perturbation, and truncates to three decimals.
`get_alpha(x, y)` separately derives the effective alpha used for reporting.

A new requirement `y >= x^(2/5)` would create additional coordination work:

- Since `y` is approximately `alpha*x^(1/3)`, the corresponding alpha floor
  is approximately `x^(1/15)`, not `x^(2/5)`.
- Three-decimal alpha truncation and the subsequent conversion to integer
  `y` can round a nominal boundary down.
- Clamping only `y` leaves a callable alpha getter that can still return a
  value below the new requested floor. A later caller could reconstruct
  an inconsistent `y` or apply the tuning perturbation twice.
- A downward double-check perturbation followed by a new lower clamp can
  select the same `y` as the original computation, reducing the intended
  independence of the check for additional inputs.
- A conservative bound such as `(iroot<5>(x)+1)^2` overshoots the required
  real bound, especially for small inputs. An exact rational-power routine
  would add arithmetic that C1 removal does not need.
- At `x <= 10^31`, both `x^2` and `y^5` near the proposed boundary can be
  about `10^62`; neither comparison operand fits in 128 bits. Their use
  in mathematical proofs does not justify evaluating them directly in C++.

Because no new floor is required, these are avoidable changes. The cleanest
C1 removal preserves the meanings of both existing alpha functions and all
seven setup paths. Existing duplicate parameter setup could be centralized
in a separate change if desired, but it is not needed to establish correctness
here.

If a mandatory `x^(2/5)` floor is later requested as an independent tuning
policy, a single `get_y_gourdon(x, requested_alpha)` boundary would be cleaner
than exposing independently callable raw-alpha selection and `y` normalization.
Raw tuning should then become private, effective-alpha reporting should use
the finalized `y`, and double-check behavior must be designed around the
final feasible range. That alternative is not the recommended C1 change.

## 9. Proposed implementation, after review and authorization

### 9.1 Minimal change

1. Delete `C1()`, `C1_libdivide()`, and `C1_arm_sve()`.
2. Delete the C1 dispatch blocks, `max_c1`, and C1-only intermediates from
   all three `AC_OpenMP_*` implementations.
3. Remove their `DBand.hpp` includes once no symbols from that header remain.
4. Update comments that list C1 as work performed by the shared thread team.
   A compact explanation that D already includes these rows is appropriate.
5. Preserve all existing AC and D fixtures and expected results, including
   low-`y` cases. Correct the contrary statements in the design document.
6. Add the required ChangeLog entry for the eventual implementation and
   perform the cleanup pass required by `AGENTS.md`.

Keep the A and C2 formulas, their row bounds, and their names. In particular,
`pi_root3_xy`, `pi_sqrty`, `xy`, `pi_sqrt_low`, the prime vectors, and the
libdivide prime vector still have uses outside C1. `sqrty` can be folded into
`pi[isqrt(y)]` if it becomes an otherwise unused intermediate.

Keep `LoadBalancerAC`, `SegmentedPiTable`, the segment range, and the
below-`y` scheduling behavior. A/C2 still have work there, including reflected
terms. In D, retain the kernel cutoff's `low >= y` condition.

This removes code and small setup/dispatch overhead. It does not eliminate
a currently executed C1 summation, so a large new runtime improvement should
not be inferred from the number of deleted lines.

### 9.2 Optional adjacent simplification in D

Section 4 proves `pi_band <= pi_root3_xy`, throughout the current valid
parameter domain. Therefore D's easy-band starting index can change from

```cpp
max3(k, pi_root3_xy, pi_band) + 1
```

to

```cpp
max(k, pi_root3_xy) + 1
```

This is valid for existing explicit-`y` D tests too, not only for a hypothetical
new clamped path. It may be included as a small adjacent simplification or
reviewed separately. If it is included, search for every remaining
`band_min_prime()` use: AC removal plus this simplification may leave only
test callers. Decide deliberately whether to remove the helper and revise
the redundant classification checks in `test/gourdon/DBand.cpp`.

### 9.3 Documentation corrections

The original design document needs these corrections after independent review:

- Section 1: C1 computes nothing for all valid current parameters.
- Section 3.4: below the threshold, there are neither band rows nor original
  C1 rows; the claimed residual C1 interval is empty.
- Sections 4.5 and 10: replace the proposed lower clamp with the stronger
  integer proof and unreachable-code removal.
- Section 10.2: the low-input `C1 empty` entries marked `no` are incorrect.
- Section 10.3: retain low-`y` fixtures; do not regenerate them merely because
  C1 is deleted. The D starting-index simplification is unconditional within
  the current valid domain.

The earlier TeX paper describes the predecessor algorithm. Updating it to
describe the full extension is broader work and should not be mixed into
this small code removal without a separate decision.

## 10. Verification performed and its limits

The analysis used independent derivations and three numerical approaches.
The counts below describe checks of bounds, not new full prime-counting runs.

1. A scratch C++14 harness called the actual production `band_min_prime()`
   and compared it with independent exact integer-root references. With
   `-O2 -Wall -Wextra -pedantic -Werror`, it passed:
   - 33,558,528 exhaustive pairs `1 <= y <= N <= 8192`;
   - 21,032,169 exhaustive pairs `2 <= x <= 100000`, `1 <= y <= isqrt(x)`;
   - 900,000 endpoint/random cases through `10^31`, respecting the quotient
     guard, with 128-bit support;
   - the same exhaustive sets and another 900,000 cases through `10^18`
     when compiled with `-DDISABLE_INT128`.
2. An independent JavaScript check passed all 21,032,170 pairs
   `1 <= x <= 100000`, `1 <= y <= isqrt(x)`. The extra case is `x=y=1`.
3. An independent C# arbitrary-precision check passed 775,721 pairs,
   combining exhaustive small inputs, 100,000 random cases through `10^31`,
   and power boundaries involving `p^4`, `p^5`, `p^6` and `y=p^2,p^3`
   with neighboring integers.

The C++ random generator used `std::mt19937_64` with seed `0x9d8ca1762`.
It generated 300,000 inputs across decimal scales and checked three `y`
choices per input: the smallest currently feasible value, `isqrt(x)`, and
a random value between them. The small feasible value also respected
`x <= y*2^62`. Root references used binary search and division-based
comparisons, avoiding overflow. These checks do not depend on Python.

The scratch C++ source and executables are in the ignored local build
directory as `build-agents-msys2-gcc/c1_bound_analysis.cpp`,
`c1_bound_analysis.exe`, and `c1_bound_analysis_64.exe`. They are local
analysis artifacts and are not part of this document's tracked changes;
the independent review must not depend on receiving them.

No implementation has been changed or built as part of the proposed C1
removal. This document does not claim that the full D identity, all compiler
targets, or all future edited code have been verified by these bound sweeps.

### 10.1 Reproduction outline independent of the scratch files

An independent checker can use this specification:

```text
exact_root(n, d):
    binary-search the largest integer r with r^d <= n
    use division comparisons, or arbitrary-precision integers

check(N, y):
    require 1 <= y <= N
    R = exact_root(N, 3)
    Q = exact_root(y, 3)
    S = exact_root(y, 2)
    reference = min(max(Q, N/y), S)
    compare production band_min_prime(N, y) with reference
    assert reference <= R

for N in 1..8192:
    for y in 1..N:
        check(N, y)

for x in 1..100000:
    for y in 1..exact_root(x, 2):
        check(x/y, y)
```

For row-level checks, independently enumerate primes `p <= isqrt(y)`:

```text
candidate_C1 = p > exact_root(x/y, 3)
has_band = exact_root(x/p, 2) < y and p*p > y/p
mixed = (x/p)/(p*p) < y

assert candidate_C1 implies has_band and mixed
assert candidate_C1 implies p <= exact_root(x, 4)
assert p <= band_min_prime(x/y, y) implies not candidate_C1
```

The `k` restriction can then be added; it only removes rows. Use exact
integers for any real-threshold cross-check. For large samples, preserve
the current quotient guard when narrowing `x/y` to `int64_t`.

## 11. Independent review and subsequent implementation checks

Claude should verify or refute each of these points against the source:

1. The three current AC implementations use the same C1 bounds shown in
   Section 4; no other C1 entry point or direct helper test executes it.
2. The inequality in Section 4 holds with floors and at equality boundaries,
   without a hidden assumption `y >= x^(2/5)`.
3. The row proof in Section 5 establishes both band eligibility and the
   mixed classification. D's upper prime cap and `k` do not omit rows.
4. The branch already includes these C1 contributions in D, with the
   correct signs, boundaries, and split between easy terms and sieve terms.
5. All PiTable arguments and prime indices in D, including the reflected
   sums and `phi_vector()` lookahead, fit the current allocations. A larger
   sieve position must not be mistaken for a larger PiTable argument.
6. The unused prefix entries discussed in Section 7.3 are not consumed
   outside the closed formula's mathematical domain.
7. All seven parameter construction paths preserve the existing bounds;
   tiny inputs, the maximum input, and the quotient guard are accounted for.
8. Existing low-`y` fixtures remain useful and need no changed expectations.
   Keep cases on both sides of the threshold and cases where the original
   C1 range would have been nonempty before band transfer.
9. The optional D simplification follows from the same theorem. If accepted,
   remove or retain the cutoff helper deliberately after checking all uses.
10. No change to alpha semantics, double-check tuning, A/C2 bounds, or
    segmentation is needed merely to remove the unreachable C1 code.

After independent agreement and user authorization, run the required checks
from `AGENTS.md` using separate Release and Debug build directories:

- Configure with tests enabled and warnings treated as errors for C and C++.
- Build, run `ctest`, and inspect any failures with `--output-on-failure`.
- Run `primecount --test` and verify `primecount 1e18 -s` returns exactly
  `24739954287740860`.
- Retain AC, D, DBand, and alpha regression coverage, and compare full
  Gourdon results with an independent implementation for selected explicit
  alphas and thread counts.
- Exercise both 64-bit and 128-bit instantiations. Cover scalar and libdivide
  AC paths; verify ARM SVE changes with the available platform or compile
  checks, reporting execution limits explicitly.

Use Debug assertions to catch table and vector range errors. The standalone
cutoff proof should be supplemented by existing direct band/reference tests,
not used as a substitute for checking the eventual edited program.
