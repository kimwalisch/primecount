# Three-prime extension of Gourdon's D hyperbola method

Specification and analysis as of commit `a8015b15` (2026-10-25). The full three-prime extension is implemented in this revision. This document covers its derivation, implementation constraints, verification, measured performance, and relationship to the C2 easy special leaves.

The original inputs were `../math_papers/Hard-Special-Leaves-Hyperbola.md` and `../math_papers/Hard_Special_Leaves_Cancellation_2.tex`, with those paths relative to the repository root. These companion documents are outside this repository and describe an earlier two-prime implementation and older default parameters. The derivation needed for the three-prime extension is included below. Related repository papers are [Hard-Special-Leaves.pdf](Hard-Special-Leaves.pdf), [Hard-Special-Leaves-SIMD-Filtering.pdf](Hard-Special-Leaves-SIMD-Filtering.pdf), and [Easy-Special-Leaves.pdf](Easy-Special-Leaves.pdf), especially Section 6 for C2.

**The extension preserves the existing cheap three-pi expression.** It adds triple-product events and a corrected boundary prefix. Independent direct-sum checks and the production test suites pass. A fixed-parameter benchmark measured a median reduction of about 0.8%, close to timing variation; the much larger reduction in query counts within the added bands does not imply a comparable reduction in total D time.

Throughout, a “three-prime leaf” means the leaf parameter m has three prime factors. The complete denominator p*m then has four prime factors. The two-prime baseline handled at most two factors in m within its hyperbola range. The extended hyperbola method handles up to three; the existing D factor-table kernel already supports more factors outside that range. Prime factor counts include multiplicity when discussing the auxiliary events, although actual Moebius leaves are squarefree.

For future implementation work, follow [AGENTS.md](../AGENTS.md), including its required paper reading, compatibility rules, preservation of existing changes, and full verification procedure. The tables and timings in this document are snapshots of the revision above; inspect the current source before changing its bounds or tuning parameters.

All quotient arguments of prefix functions are rounded down, including those where floors are omitted for readability. Write P^-(n) for the least prime factor, with P^-(1)=infinity, and Omega(n) for the number of prime factors counted with multiplicity. The supported input range is x<=10^31, with y^2<=x. D's caller invokes the hyperbola helper only for b>k and p^2<=y; extending its internal guard does not remove those caller assumptions.

For a fixed sieving prime p=p_b, write

$$
X=\lfloor x/p\rfloor,\quad H=\lfloor\sqrt X\rfloor,\quad
L=\lfloor X/y\rfloor,\quad \Phi(t)=\phi(t,b-1),
$$

$$
a(n)=[P^-(n)\ge p],\qquad c(n)=\mu(n)[P^-(n)>p].
$$

Thus Phi(t)=sum_{n<=t} a(n) and M(t)=sum_{n<=t} c(n). The distinction between factors >=p in a and factors >p in c is essential for repeated products and products containing p.

The contribution to be rearranged is

$$
D_p=-\sum_{L<m\le y}c(m)\Phi(\lfloor X/m\rfloor).
$$

The two-prime baseline handles H<y<p^3. The extension adds

$$
\boxed{H<y,\qquad p^3\le y<p^4.}
$$

As usual, y^2<=x. In these new cases, L>=p^2: x>=y^2 and y>=p^3 imply X>=p^2*y. Also x>=p^3*y, so all leaves m<=y of these primes are hard. This extension does not require moving further easy leaves out of AC.

**The weights in the events vector are cancellation coefficients.** They are a(n)+c(n), rather than the number of prime factors. The distinction becomes particularly useful for triples:

| Triple n, all factors >=p | a+c, usual reflected-query coefficient | a-c, upper correction coefficient |
| --- | ---: | ---: |
| Three distinct primes, all >p | 0 | 2 |
| A repeated prime, all factors >p | 1 | 1 |
| At least one factor equal to p | 1 | 1 |

Thus no event weight 3 is needed. Generic squarefree triples cancel from the usual reflected sieve-query sum, just as primes do. They still contribute cheap terms and, when above H, positive correction queries.

**Simply changing p^3 to p^4 in the eligibility guard would be incorrect.** Retain the expression

$$
Q(t)=2\pi(t)+\pi(\lfloor t/p\rfloor)+\pi(\lfloor\sqrt t\rfloor)-4b+2.
$$

For t>=p^2 it is exactly the contribution of integers with at most two prime factors to Phi(t)-M(t). It ceases to be the entire difference at p^3.

Let Omega count prime factors with multiplicity. For t<p^4 define

$$
T_3(t)=\#\{qrs\le t:p\le q\le r\le s,\ q,r,s\text{ prime}\},
$$

and let S_3(t) count triples with p<q<r<s. Let M_2(t) be the restricted Mobius prefix with at most two factors, which is exactly what the existing descending recurrence computes:

$$
M_2(t)=1-\pi(t)+b+
\sum_{b<i\le\pi(\sqrt t)}(\pi(\lfloor t/p_i\rfloor)-i).
$$

Then, for p^2<=t<p^4,

$$
M(t)=M_2(t)-S_3(t),\quad
\Phi(t)-M(t)=Q(t)+T_3(t)+S_3(t).
$$

The last equality explains the obstacle to a naive extension: its generic triple layer survives with coefficient 2. Evaluating that layer separately at every argument X/n would introduce expensive nested prime sums. That cost can be avoided by transposing the entire correction.

**An exact three-prime formula keeps Q cheap.** Define the auxiliary prefix

$$
\widetilde M(t)=\Phi(t)-Q(t)=M_2(t)+T_3(t),
$$

which is not the restricted Mobius prefix once triples occur. Define d_3(k)=a(k)-c(k) on the triple products; its values are 1 or 2 as in the table. The new formula is

$$
\boxed{\begin{aligned}
D_p={}&[\Phi(H)-Q(H)]\Phi(H)
       -[M_2(y)+T_3(y)]\Phi(L)\\
&+\sum_{L<n\le H}a(n)Q(\lfloor X/n\rfloor)\\
&-\sum_{L<n\le H}[a(n)+c(n)]\Phi(\lfloor X/n\rfloor)\\
&+\sum_{\substack{H<k\le y\\\Omega(k)=3\\P^-(k)\ge p}}
       d_3(k)\Phi(\lfloor X/k\rfloor).
\end{aligned}}
$$

The n sums now include all eligible triples below H, as well as primes and semiprimes. Each such triple adds Q(X/n). A triple with repeated factors or a factor p also needs a negative weight-1 query; a squarefree triple entirely above p needs no query there. Triples above H instead need the final positive correction query.

The boundary change from M_2(y) to M_2(y)+T_3(y) is essential. When y>=p^3, M_2(y) is not the true restricted Mobius prefix; the old guard excluded those entries from use. Lowering the guard without changing this boundary would give incorrect answers.

A short derivation uses an arbitrary-weight hyperbola identity. For prefix sums A(t)=sum_{n<=t} a(n) and C(t)=sum_{n<=t} c(n), splitting at H and reversing the part with m>H gives

$$
-\sum_{L<m\le y}c(m)A(X/m)
=C(H)A(H)-C(y)A(L)
-\sum_{L<n\le H}[c(n)A(X/n)+a(n)C(X/n)].
$$

This identity depends only on the prefix sums and H<y, with L=floor(X/y); it places no arithmetic restriction on the weights. Define e_2(n)=(a-c)(n) for Omega(n)<=2 and zero otherwise, so its prefix is Q at arguments >=p^2. Put c_hat=a-e_2. On integers with at most two factors, c_hat=c; on triples with a=1, c_hat=1. Thus its prefix is M_tilde=Phi-Q. Apply the identity to a and c_hat, for which the cheap difference is Q. The desired sum differs from that auxiliary sum by

$$
\sum_{\substack{L<k\le y\\\Omega(k)=3}}d_3(k)\Phi(X/k).
$$

For k<=H this combines with the auxiliary coefficient -(a+c_hat)=-2 to give the desired coefficient -(a+c). For k>H it remains as the final positive correction sum. This proves the boxed identity without evaluating a triple-count prefix at each n.

**The current sorter representation remains sufficient.** The query positions have this order:

| Work | Query positions | Query coefficient in the sum |
| --- | --- | --- |
| Boundary at L | L | -(M_2(y)+T_3(y)) |
| Triple corrections k>H, decreasing k | [L,H] | +1 or +2 |
| Boundary at H | H | Boundary product |
| Semiprimes and noncancelling triples n<=H, decreasing n | >=H and <y | -1 or -2 |
| Remaining factor-table leaves m<=L | >=y | Existing Mobius coefficient |

Positive and negative batches can use separate calls to the existing sorter; signs need not be encoded as new weights. Unique sorted prime factorizations q<=r<=s preserve the sorter's distinct-product invariant, and products remain coprime to 30 when p>5. Triple cheap terms with zero query weight need no sorting. Products sharing the same integer quotient X/k could optionally be combined, but the censuses below do not assume that optimization.

Generate triples directly from ordered primes, restricting the last prime using integer quotient bounds for the current segment. Do not regenerate every triple up to y on every segment. The positive batch uses H<k<=y and positions X/k in the segment; the negative batch uses L<n<=H with positions X/n in the segment. Keep the current max_m=L lowering and the remaining factor-table work. The existing native SIMD dispatch can remain in place.

Precompute T_3(y) once per prime, alongside the existing M_2 recurrence. If p_next is the next prime after p, count by the least factor:

$$
T_{3,p}(y)-T_{3,p_{next}}(y)
=\sum_{p\le r\le\sqrt{y/p}}
  [\pi(\lfloor y/(pr)\rfloor)-\pi(r)+1].
$$

Descending p therefore needs only semiprime sums for the added prefix. It does not require enumerating every triple just to compute a boundary value. Retain the current unsigned accumulation and use quotient comparisons for eligibility and product bounds, rather than forming overflowing p^4 values in generic code. Triple products that are actually generated are <=y and fit the existing 64-bit event field.

**The production implementation maps these formulas to the following files.**

| File | Responsibility |
| --- | --- |
| [DHyperbola.cpp](../src/gourdon/DHyperbola.cpp) | Eligibility, boundary terms, positive triple corrections, ordinary semiprime events, and lower triple terms |
| [DHyperbola.hpp](../src/gourdon/DHyperbola.hpp) | Shared prime-sum helper and old easy terms with query positions below p^2 |
| [D.cpp](../src/gourdon/D.cpp) | Precompute M_2(y)+T_3(y) for eligible primes, thread scheduling, and dispatch |
| [SemiprimeEventSorter.hpp](../src/gourdon/SemiprimeEventSorter.hpp) | Sort distinct products with weights 1 or 2; its name is retained although some products are triples |
| [test/gourdon/DHyperbola.cpp](../test/gourdon/DHyperbola.cpp) | Independent dense references for each segment, query order, boundary placement, and integer widths |
| [test/gourdon/D.cpp](../test/gourdon/D.cpp) | Whole-D reference comparisons and transition cases |

The five `D_*.hpp` kernels retain their existing helper interface. The triple prefix only needs entries with

$$
\max\bigl(\lfloor y^{1/4}\rfloor,\lfloor x/y^2\rfloor\bigr)<p
\le\lfloor y^{1/3}\rfloor,\qquad b>k.
$$

Compute floor(x/y^2) as `xy / y`, where `xy = x / y`. Descend from pi(floor(cuberoot(y))) with the triple count initially zero. For each p, add the least-factor-p contribution to the triple count **before** adding that count to the stored M_2(y) entry. M_2 excludes the prime p from its factors, while T_3 includes it; confusing those conventions changes the boundary value.

In `hyperbola_sum()`, reject H>=y before doing the extra eligibility divisions. After `square = p*p`, the safe comparisons `square <= y/p` and `square <= y/square` detect y>=p^3 and y>=p^4, respectively. Reject the latter. Do not replace these checks by forming p^4 in a 64-bit temporary.

For a sieve segment [low,high), use `query_low = max(low,p*p)`. Positive correction products lie in

$$
\max(H,\lfloor X/high\rfloor)<n
\le\min(y,\lfloor X/query\_low\rfloor).
$$

The lower events and cheap terms lie in

$$
\max(L,\lfloor X/high\rfloor)<n
\le\min(H,\lfloor X/query\_low\rfloor).
$$

Clamp wide quotients before narrowing them to 64 bits, using the project's mixed-width min helpers. Enumerate q=primes[i] up to cuberoot(max_n), then r=primes[j] with j>=i up to sqrt(max_n/q). For each pair qr, the last-prime indices run from max(j,pi(floor(min_n/qr))+1) to pi(floor(max_n/qr)). These bounds ensure every generated product is in the current range and fits its 64-bit destination.

For upper triples, enumerate all allowed last-prime indices k and use weight 2 precisely when i>b and i<j<k; otherwise use weight 1. Sort these events and add their queries before querying H.

For lower triples, products with q=p or q=r generate weight-1 events for every allowed last prime. If q>p and q<r, emit the r=s event only when it is in range. For the generic triples q<r<s, set

$$
\alpha=\max(\lfloor min\_n/(qr)\rfloor,r),\qquad
\beta=\lfloor max\_n/(qr)\rfloor.
$$

When alpha<beta, reuse `sum_phi_minus_mertens(X/(qr), p, b, alpha, beta, ...)` to sum Q over the last prime. Its inversion identities remain valid when the arguments exceed p^3: they sum the expression Q, not the full difference Phi-M. Keep generic triples out of the event vector. The existing final event loop supplies Q once for each semiprime or degenerate triple event; adding Q separately for those triples would double-count it.

The new primes satisfy L>=p^2, so their contributions need no change to `easy_hyperbola_sum()`. That function still handles the old easy cases with y<p^3. At x<=10^31, y<=sqrt(x)<3.163*10^15. Generated products and individual prefixes fit signed 64 bits, but X and accumulated products may need 128 bits. Preserve the existing unsigned modular accumulation and the C++14 implementation/C++11 public-header requirements.

**The restricted prototype provides a useful intermediate implementation.** Use only

$$
\boxed{p^7>x,\qquad p^3\le y.}
$$

Given y^2<=x, these imply H<p^3, H<y<p^4, and L>=p^2. There are no triple n below H. The existing prime and semiprime loops can therefore remain as they are; only the adjusted boundary and positive triple events are added. This restricted version was implemented and benchmarked before the full version with triples on both sides of H. It covers 13 new primes at x=10^20; the full extension covers 75.

**Use the branch's current alpha, rather than the documents' old default tables.** Geometrically, new primes must satisfy approximately

$$
\max(x/y^2,y^{1/4})<p\le y^{1/3}.
$$

In particular y>x^(3/7) is necessary. If y<=x^(3/7), any prime with H<y already satisfies y<p^3, so this extension opens no new ranges. For the current source's default tuning, the necessary condition does hold at the following representative inputs:

| x | Current default y | Newly eligible primes | Restricted H<p^3 primes |
| --- | ---: | ---: | ---: |
| 10^14 | 1,224,474 | 9 | 2 |
| 10^16 | 8,720,917 | 14 | 2 |
| 10^18 | 67,144,000 | 32 | 6 |
| 10^20 | 551,077,176 | 75 | 13 |
| 10^24 | 34,605,200,000 | 315 | 72 |
| 10^31 | 30,157,949,145,395 | 2,017 | 410 |

For example, at 10^20 the documents list y=245,535,363, whereas the current alpha coefficients produce y=551,077,176. With the old value there would be no new range. With the current value there are 75 added primes. Exact integer guards, rather than these approximate root bounds, should control an implementation.

**The workload census measures query reduction, not execution time.** These counts cover only the newly eligible bands L<m<=y, comparing the two-prime baseline's factor-table leaf queries with the three-factor formula. They count one query per event product and two boundaries per prime, without deduplicating equal positions.

| x | Two-prime baseline queries in the added bands | Three-prime formula queries | Reduction within these bands |
| --- | ---: | ---: | ---: |
| 10^14 | 261,156 | 50,135 | 80.8% |
| 10^16 | 2,788,987 | 533,873 | 80.9% |
| 10^18 | 52,573,982 | 9,862,923 | 81.2% |
| 10^20 | 1,351,576,213 | 242,882,561 | 82.0% |

The 10^20 replacement comprises 240,371,876 semiprime queries, 149,656 noncancelling lower triple queries, 2,360,879 upper triple correction queries, and 150 boundary queries. The lower triples also contribute 1,512,832 mathematical Q terms, including those that need no query; the implementation aggregates generic terms through inversion instead of evaluating all of them individually. The usual prime and semiprime Q work remains and is not included in the query-count column.

At 10^20, the restricted H<p^3 prototype alone changes 345,016,213 original queries into 54,152,686 queries. It needs only 11,155 upper triple correction queries, in addition to its existing-form semiprime queries and boundaries.

These figures are not predictions of an 81-82% speedup in D. They apply to additional bands only, and account for neither pi lookups nor generation, sorting, segmentation, sieving, load balancing, or cache behavior. Benchmark D at a fixed y first; only after that evaluate total prime-counting time and possible retuning.

**Measured performance is modest.** The two-prime baseline at `e68f53f3` and the full extension were built with Homebrew LLVM/Clang 23.1.1, Release optimization, OpenMP, and `-Wall -Wextra -pedantic -Werror` for both C and C++. D selected ARM NEON on the test machine. Both binaries used

```sh
./primecount 1e19 -D --time --alpha=89.048
```

This fixes the default-at-that-revision y=191848038. Every run returned D=197179281410717181. One warmup per binary preceded six measured pairs with alternating baseline/candidate order. Builds and other tests did not run during the timed comparisons. The reported wall times include initial PiTable work.

| Pair | Run order | Baseline seconds | Full extension seconds |
| --- | --- | ---: | ---: |
| 1 | Baseline, extension | 9.521 | 9.600 |
| 2 | Extension, baseline | 9.816 | 9.673 |
| 3 | Baseline, extension | 9.933 | 9.974 |
| 4 | Extension, baseline | 10.350 | 10.001 |
| 5 | Baseline, extension | 10.205 | 10.002 |
| 6 | Extension, baseline | 10.226 | 10.104 |
| Median | | 10.069 | 9.988 |
| Mean | | 10.008 | 9.892 |

The median reduction is about 0.8%, and the mean reduction about 1.2%. Four of six pairs favor the extension, but the difference is close to timing variation. These measurements establish neither a large improvement nor a general speedup across machines or input sizes.

The earlier restricted prototype had no reliable measured gain: its phase's baseline/candidate medians were 9.922/10.068 seconds and means 10.370/9.984 seconds. That phase included a large baseline outlier and timing drift. It also preceded the correction that preserves the H>=y short circuit before computing the extra eligibility divisions. Compare binaries within the same alternating phase, rather than comparing absolute times across phases.

**The standalone mathematical validation is independent of the production event generator.** A C++14 checker compares the formulas against dense least-prime-factor, Mobius, and partial-sieve prefixes. It passes 538,452 checks of E=Q+the triple correction, 3,943 triple-weight checks, and 9,000 physical bands with y^2<=x and p^3<=y<p^4. Of those bands, 804 have H<p^3 and 8,196 have H>=p^3. It also checks the transposed correction and the boundary adjustment relative to the actual current recurrence. A separate independent checker passes 59,429 small cases, including 33,114 with triple indices below H.

**Production verification is recorded separately from the standalone mathematics.** At the implementation revision, both Release and Debug passed all 69 CTest tests, `primecount --test`, and `primecount 1e18 -s`, with the exact result 24739954287740860. Both used Homebrew Clang and OpenMP with warnings treated as errors. Debug also enabled assertions, `-O1 -fsanitize=undefined -fno-sanitize-recover=all`, and produced no sanitizer diagnostics.

The committed segment test preserves the original 25 cases and adds 141 accepted cases and three rejection cases. It checks both H<p^3 and H>=p^3; products containing p; repeated factors; squarefree triples; y near p^3 and p^4; H near p^3; and placement across segment sizes 240 and 1680. It exercises 64-bit and 128-bit accumulator types and 32-bit and 64-bit prime vectors. The whole-D test adds nine independent reference comparisons around triple-product transitions, using both accumulator widths.

For subsequent changes, run the complete verification required by [AGENTS.md](../AGENTS.md) and the Debug checks in [test/README.md](../test/README.md). Keep the exact rejections at y=p^4 and H>=y. Dense tests that instantiate uint128_t with small values do not prove coverage of X>2^64; when changing quotient handling, independently check wide-quotient bounds and add suitable real wide-input coverage. The mathematical bounds above complement rather than replace runtime verification.

The following historical scratch artifacts were used during the analysis. They are outside the repository, may not exist on another machine, and are not dependencies of this specification or of the committed tests:

- `/private/tmp/verify_hyperbola_three.cpp` and `/private/tmp/verify_hyperbola_three.txt`
- `/private/tmp/check-hyperbola-three.cpp`
- `/private/tmp/dhyperbola-census.cpp` and its `-1e14.txt`, `-1e16.txt`, `-1e18.txt` outputs
- `/private/tmp/census_three_1e20.cpp` and `/private/tmp/census_three_1e20.txt`
- `/private/tmp/dhyperbola-geometry.js` and `/private/tmp/dhyperbola-geometry.txt`

The dense verifier and the 10^20 census were also rebuilt with C++14, `-Wall -Wextra -pedantic -Werror`, and UndefinedBehaviorSanitizer with recovery disabled. Both completed successfully with no diagnostics and the same results.

The standalone checks establish the identities and workload counts; the committed tests and implementation runs cover production integration. Runtime improvement remains limited to the small measured change above. Alpha tuning was not changed for this extension.

**Applying the same cancellation to C2 reproduces its existing inversion.** C2's sieving primes satisfy p>sqrt(y), so every admissible cofactor m<=y is prime. The three-prime extension adds primes with p^3<=y, which are entirely below this C2 range. There is no additional three-factor layer to exploit in C2.

For one C2 prime, define

$$
R=\max(p,\lfloor X/p^2\rfloor).
$$

Its full contribution, before folding, is

$$
C_{2,p}=\sum_{\substack{R<q\le y\\q\ \mathrm{prime}}}
  [\pi(\lfloor X/q\rfloor)-b+2].
$$

The apparent upper bound min(y,floor(X/p)) in the kernel reduces to y: Gourdon's setup enforces p<=x_star<=sqrt(x/y), hence p^2*y<=x and L>=p. Also y<p^2 gives R<=L. If the row is nonempty, some q<=y satisfies X<p^2*q<=p^2*y<p^4, so H<p^2 even in rows with H>=y.

When H>=y, there is no nonempty top band with upper end y for this reflection. When H<y, the useful bounds are

$$
R\le L\le H<y<p^2.
$$

Every reflected rough integer is prime. For all the band query positions p<=t<y<p^2,

$$
\Phi(t)=\pi(t)-b+2,\qquad M(t)=1-\pi(t)+b,
$$

$$
E(t)=\Phi(t)-M(t)=2\pi(t)-2b+1=2\Phi(t)-3.
$$

The coefficient a+c vanishes on reflected primes, but their remaining E term contains exactly the same pi lookup with coefficient 2. In hard D, replacing a sieve query by pi lookups changes the expensive operation. In C2, Phi already comes from a pi lookup, so this symbolic cancellation does not remove that operation.

Substitution into the weighted identity gives the exact top-band formula

$$
\begin{aligned}
\sum_{\substack{L<q\le y\\q\ \mathrm{prime}}}
  [\pi(X/q)-b+2]
={}&\pi(y)\pi(L)-\pi(H)^2\\
&-(b-2)[\pi(y)-\pi(L)]
 +2\sum_{\substack{L<q\le H\\q\ \mathrm{prime}}}\pi(X/q).
\end{aligned}
$$

This can also be proved directly by counting prime pairs qr<=X with both factors <=y. Split at H: the pair count is twice the part with q<=H, minus pi(H)^2. For q<=L, every prime r<=y is allowed, giving a rectangle pi(L)*pi(y). Subtract that rectangle to restrict q to (L,y], then subtract the constant b-2 once per prime in that interval. All floors are already accounted for by L and H.

Adding the untouched sparse range R<q<=L gives

$$
\begin{aligned}
C_{2,p}={}&\pi(y)\pi(L)-\pi(H)^2
 -(b-2)[\pi(y)-\pi(R)]\\
&+\sum_{\substack{R<q\le L\\q\ \mathrm{prime}}}\pi(X/q)
 +2\sum_{\substack{L<q\le H\\q\ \mathrm{prime}}}\pi(X/q).
\end{aligned}
$$

Thus there is one leaf-loop lookup per prime R<q<=H, with a doubled coefficient on the reflected part. This is already the algorithm in `C2()` in [AC_default.hpp](../src/gourdon/AC_default.hpp), `C2_libdivide()` in [AC_libdivide.hpp](../src/gourdon/AC_libdivide.hpp), and `C2_arm_sve()` in [AC_arm_sve.hpp](../src/gourdon/AC_arm_sve.hpp). Each has a boundary correction and three sparse-loop ranges, using multipliers 1, 2, and 1. Section 6 of [Easy-Special-Leaves.pdf](Easy-Special-Leaves.pdf) derives this additional application of Gourdon's inversion equality and the overlap with the sparse leaves.

The source uses `max_clustered_global = prev_prime(y)` and integer reflected endpoints such as floor(X/(H+1)). These adjust the boundaries and which sparse terms get multiplier 2. They do not change the fact that the three loops together visit the ordinary sparse primes R<q<=H once. Preserve this endpoint handling when comparing implementations; replacing it casually by real-valued square-root bounds can introduce diagonal errors.

A source-bound census at the same default tuning gives:

| x | Unfolded C2 leaf lookups | Current folded C2 leaf lookups | Weighted hyperbola plus remaining leaves |
| --- | ---: | ---: | ---: |
| 10^12 | 1,069,044 | 337,540 | 337,540 |
| 10^14 | 24,186,914 | 4,700,998 | 4,700,998 |
| 10^16 | 464,561,629 | 74,679,189 | 74,679,189 |

These counts are a mathematical census of the source loop bounds, including the actual last-prime endpoint. They exclude constant boundary work and segment setup; they are not production instrumentation or timings. For a folded row, count pi(H)-pi(R); for an unfolded row with H>=y, count pi(y)-pi(R). The weighted formulation has pi(L)-pi(R) remaining lookups plus pi(H)-pi(L) reflected lookups, the same total. The band and full-row identities were also checked against 6,000 direct sums using primes 7 through 199, integer floors, y^2<=x, and the C2 support bounds.

The C2 census scratch files were `/private/tmp/c2-hyperbola-census.js` and `/private/tmp/c2-hyperbola-census.txt`; they are optional historical artifacts, like the D scratch files above.

Moving the upper C2 band into D could change scheduling, cache locality, or table access. It would not supply an additional mathematical reduction from this cancellation mechanism, and sparse arguments at or above y cannot all use D's pi table bounded by y. In particular, do not extend the D caller past sqrt(y) and keep using its p^2<=y assumptions or Q formula at t<p^2. Any such experiment needs the separate prime-only formulas above and an independent performance comparison. Further mathematical savings would require another identity, rather than another application of this same per-prime reflection.
