# Three-prime hyperbola method for the hard special leaves

The hyperbola method in `src/gourdon/DHyperbola.cpp` computes the special
leaves `x/(p*y) < m <= y` of a sieving prime `p`. The original method
requires `sqrt(x/p) < y < p^3`, so these `m` contain at most two prime
factors. The extension permits `sqrt(x/p) < y < p^4`, including products
of three prime factors. Prime factors are counted with multiplicity.

The important implementation choice is to reflect the three-prime
correction as a whole. Evaluating a three-prime prefix separately for
every cancelled prime leaf would introduce expensive nested prime sums.

## Prefix functions

Fix `p = primes[b]` and put

$$
X=\lfloor x/p\rfloor,\qquad
L=\lfloor X/y\rfloor,\qquad H=\lfloor\sqrt X\rfloor.
$$

Let `lpf(1) = infinity`, and define

$$
a(n)=[\operatorname{lpf}(n)\ge p],\qquad
c(n)=\mu(n)[\operatorname{lpf}(n)>p],
$$

with prefixes $\Phi(t)=\sum_{n\le t}a(n)=\phi(t,b-1)$ and
$M(t)=\sum_{n\le t}c(n)$. Let $\Phi_2$ and $M_2$ restrict the same
prefixes to integers with at most two prime factors. For $t\ge p^2$,

$$
E_2(t)=\Phi_2(t)-M_2(t)
=2\pi(t)+\pi(\lfloor t/p\rfloor)
+\pi(\lfloor\sqrt t\rfloor)-4b+2.
$$

This expression is implemented by `phi_minus_mertens()`. It equals
$\Phi(t)-M(t)$ only when $t<p^3$.

For $t<p^4$, let $Q(t)$ count all products of three primes $\ge p$,
and let $T(t)$ count products of three distinct primes $>p$. Then

$$
\Phi(t)=\Phi_2(t)+Q(t),\qquad M(t)=M_2(t)-T(t),
$$

so the missing correction is

$$
C_3(t)=\Phi(t)-M(t)-E_2(t)=Q(t)+T(t).
$$

Equivalently, $C_3$ is the prefix over three-prime products $d=qrs$,
$p\le q\le r\le s$, with weight $w(d)=a(d)-c(d)$: weight 1 if
$q=p$, $q=r$ or $r=s$, and weight 2 otherwise.

## Reflecting the correction

The derivation below concerns the new primes $p^3\le y<p^4$. They have
$L\ge p^2$, as shown below, so every argument of $E_2$ is in its domain.

The weighted Dirichlet hyperbola identity gives the contribution of
the leaves $L<m\le y$ as

$$
\begin{aligned}
D_p={}&M(H)\Phi(H)-M(y)\Phi(L)\\
&+\sum_{L<n\le H}a(n)(\Phi(X/n)-M(X/n))\\
&-\sum_{L<n\le H}(a(n)+c(n))\Phi(X/n).
\end{aligned}
$$

All function arguments are rounded down. Substitute
$\Phi-M=E_2+C_3$. Interchanging the sums in the additional term gives

$$
\begin{aligned}
\sum_{L<n\le H}a(n)C_3(X/n)
={}&C_3(H)\Phi(H)-C_3(y)\Phi(L)\\
&+\sum_{\substack{H<d\le y\\\Omega(d)=3}}
w(d)\Phi(X/d).
\end{aligned}
$$

Indeed, a triple $d\le H$ permits every $n\le H$, since $H^2\le X$.
A triple $H<d\le y$ permits precisely $n\le\lfloor X/d\rfloor\le H$.
The lower endpoint contributes $\Phi(L)$ in both cases, since
$\lfloor X/d\rfloor\ge L$ for $d\le y$.

Since $M(t)+C_3(t)=\Phi(t)-E_2(t)=M_2(t)+Q(t)$, the final formula is

$$
\begin{aligned}
D_p={}&(\Phi(H)-E_2(H))\Phi(H)
-(M_2(y)+Q(y))\Phi(L)\\
&+\sum_{L<n\le H}a(n)E_2(X/n)\\
&-\sum_{L<n\le H}(a(n)+c(n))\Phi(X/n)\\
&+\sum_{\substack{H<d\le y\\\Omega(d)=3}}
w(d)\Phi(X/d).
\end{aligned}
$$

For primes $n>p$, the sieve coefficient $a(n)+c(n)$ vanishes. It also
vanishes for products of three distinct primes $>p$. A three-prime
product containing $p$ or a repeated factor has sieve coefficient 1.
The semiprime coefficients remain 1 for $pq$ and $q^2$, and 2 for
distinct primes $q,r>p$.

## Boundary coefficients and implementation

`generate_hyperbola_coefficients()` precomputes $M_{2,b}(y)+Q_b(y)$.
The existing two-prime recurrence, starting at $R=\pi(\sqrt y)$, is

$$
M_{2,R}(y)=1-\pi(y)+R,\qquad
M_{2,b-1}(y)=M_{2,b}(y)+\pi(\lfloor y/p_b\rfloor)-b-1.
$$

For the triples, descending through the primes adds products whose
least prime factor is exactly $p_b$:

$$
Q_b(y)=Q_{b+1}(y)
+\sum_{i=b}^{\pi(\sqrt{y/p_b})}
\left(\pi\!\left(\left\lfloor\frac{y}{p_b p_i}\right\rfloor\right)-i+1\right).
$$

Start with $Q_{\pi(\sqrt[3]{y})+1}(y)=0$. The inner sum counts the
remaining factors $p_i\le p_j$, including repeated primes. Only primes
eligible for the three-prime hyperbola method need these coefficients.

`triple_events()` generates triples in their unique order
$p\le q\le r\le s$. Above $H$ it generates the reflected sieve queries
with weights 1 or 2. Below $H$ it appends the triples with sieve
coefficient 1 to the semiprime events. For the other triples it sums
$E_2(X/(qrs))$ directly over the third prime, with independent divisions
batched four at a time. Reflecting these small third primes would expand
their interval, so inversion is reserved for the main prime terms.

Within each sieve segment, queries occur in this order:

1. The boundary at $L$.
2. Reflected triples $d>H$, sorted by decreasing $d$.
3. The boundary at $H$.
4. Semiprimes and the remaining triple events $n\le H$, sorted by
   decreasing $n$.
5. The ordinary leaves $m\le L$, processed by the factor table.

The corresponding positions increase from $L$ through $H$ to positions
below $y$, followed by positions at least $y$. Each term is evaluated
only in the segment containing its position. Equality of consecutive
positions is allowed.

## Bounds and the easy leaves

The new primes satisfy $p^3\le y<p^4$. Since $y^2\le x$, they satisfy
$L\ge\lfloor y/p\rfloor\ge p^2$ and $x/p^3\ge y$. Consequently all
their leaves are hard, and the C1 hyperbola terms in AC are unchanged.
For the existing primes with $y<p^3$, $Q(y)=0$ and the original formulas
apply, including the separate formulas for positions below $p^2$.

For the supported range $x\le10^{31}$, $y\le\sqrt x<3.17\cdot10^{15}$.
All triple products, prefix counts and $E_2$ values fit signed 64-bit
integers. Products are formed only after division bounds establish
$qrs\le y$. The test $y<p^4$ uses `p*p > y/(p*p)`, without forming
$p^4$. All table arguments and hyperbola positions are at most $y$.
The query lower bound is at least $L$. Thus divisions whose results
are narrowed to 64 bits are bounded by $X/L<2y$ when $L>0$; for $L=0$
we have $X<y$. Merely dividing $X$ by $p^2$ before clamping could
overflow 64 bits in the new three-prime range.

The modified boundary prefix $M_2(t)+Q(t)$ has weights of magnitude at
most 1 on disjoint two-prime and three-prime sets, so its absolute value
is at most $t$. Hence each boundary product has magnitude at most
$yL\le X$ or $H^2\le X$. The existing unsigned 64-bit or 128-bit
accumulation handles intermediate sums modularly; the final signed
result is recovered as before.

## Verification and benchmarks

On the dual-socket AMD EPYC 7642 server (96 cores, 192 threads), using
GCC 16.1.0, Release mode and `-mpopcnt`, the isolated sequential runs
on 2026-10-31 gave:

| Command | Original | Three-prime extension | Time reduction |
| --- | ---: | ---: | ---: |
| `primecount 1e22 -D --time` | 60.119 s | 57.834 s | 3.8% |
| `primecount 1e23 -D --time` | 230.076 s | 220.677 s | 4.1% |

Both builds returned `160597612714845461240` at `1e22` and
`1532523717645331090539` at `1e23`. The compiler flags also included
`-Wall -Wextra -pedantic -Werror`, with warnings treated as errors for C
as well. No builds or correctness tests overlapped these timing runs.

Release and Debug builds passed all 68 CTest tests and `primecount
--test`. The Release build returned `24739954287740860` for `primecount
1e18 -s`. An additional 32 comparisons against the original D build
passed with different tuning factors and thread counts.

`test/gourdon/DHyperbola.cpp` checks each segment against independent
dense partial-sieve and Moebius prefixes, and checks the total against
the original leaves. It covers repeated and distinct prime factors,
cube and fourth-power boundaries, segment boundaries, coefficient
cutoffs, both prime-vector widths, and 64-bit and 128-bit arithmetic.
