Computing special leaves
========================

This file describes how to implement the algorithm to find the
special leaves (sub-algorithm of the Lagarias-Miller-Odlyzko prime
counting algorithm) as described in section 5 (pages 556-557) in the
paper "Computing pi(x) The Meissel-Lehmer Method", Mathematics of
Computation, 44 (1985), pp. 537–560, by J. C. Lagarias, V. S. Miller
and A. M. Odlyzko.

Definitions and bounds
----------------------

* pi(x) counts the primes below x.
* lpf(n) denotes the smallest prime factor of n.
* mu(n) is the Möbius function
* x^(2/5) > y > x^(1/3)
* j = pi(sqrt(y))
* a = lower bound of the current segment (algorithm uses the segmented sieve of Eratosthenes)
* b = upper bound of the current segment (...)
* k is a small constant e.g. 5, used to pre-sieve the multiples of the first k primes

Prerequisities
--------------

* Calculate for each prime pk = sqrt(y) two parallel tables Ak and Mk.
* The value of Ak(j) is the jth square-free n = y such that lpf(n) = pk.
* The value of Mk(j) is mu(Ak(j)).
* For each prime pk = sqrt(y) we store a table called Nk such that l = Nk(j) satisfies Ak(l - 1) < j*pk = Ak(l).

<p>These tables may be computed in time O(y log x) and take up
space O(y log log x).</p>

### Notes
* _**parallel table**_ means two dimensional data structure e.g. A[k][j].
* Mk(j) does not seem to be used in the algorithm?!

Ak(j)
-----

Calculate for each prime pk = sqrt(y) the parallel table Ak. The value
of Ak(j) is the jth square-free n = y such that lpf(n) = pk. Ak(j) can
be calculated efficiently using the code below:

```C++
// A[k][j]
vector<vector<int> > A(/* size = */ pi[sqrt(y)] + 1, vector<int>(1));

// the square_free[n] table returns true if n is square free
// the pi[n] table counts the primes below n
// the lpf[n] table returns the least prime factor of n
// these 3 tables must all be of size y

for (int n = 1; n <= y; n++)
    if (square_free[n] && lpf[n] <= sqrt(y))
        A[pi[lpf[n]]].push_back(n);
```

### Notes

In the <a href="#algorithm">algorithm for finding the special leaves</a>
```if (j >= A[k].size())``` (line 4) then the jth square-free n > y,
in this case we can continue with the next loop iteration.

Nk(j)
-----

For each prime pk = sqrt(y) we store a table called Nk such that
l = Nk(j) satisfies Ak(l-1) < j*pk = Ak(l). Our algorithm simply
iterates over all indices of the Ak table and for each Ak(l) finds
a j which satisfies Ak(l-1)/pk < j = Ak(l)/pk. We use a vector of maps
as our data structure, querying N(k, j) uses O(log n) operations.

```C++
// N[k][j]
vector<map<int, int> > N_maps(A.size());

for (int k = 1; k < A.size(); k++) {
    for (int l = 1; l < A[k].size(); l++) {
        int j = A[k][l] / primes[k];
        N_maps[k][j] = l;
    }
}

int N(int k, int j) {
    int l = N_maps[k].lower_bound(j)->second;
    return l;
}
```

### Notes

In the <a href="#algorithm">algorithm for finding the special leaves</a>
we have to add bounds checking in line 5.

Algorithm
---------

<p>To find all special n ? [a, b) we use these two tables together with
the procedure from page 557 (figure 2).</p>

```C++
// Special leaves with lpf(n) <= sqrt(y)
for i := k to pi(sqrt(y))
    for j := i + 1 to pi(sqrt(y))
        l :=  Nj((a - 1) / primes[j] + 1)
        while Aj(l) <= b / primes[i]
            // it is a special leaf
            process(primes[i] * Aj(l));
            l := l + 1

// Special leaves which are the product of two primes
for i := pi(sqrt(y)) to pi(y)
    l := pi(a / primes[i]) + 1;
    while primes[l] <= b / primes[i]
        // it is a special leaf
        process(primes[i] * primes[l]);
        l := l + 1
```

### Notes

In the first revision of the paper "Computing pi(x) The Meissel-Lehmer
Method" the first part of the above algorithm contains an error, Aj(l)
and Nj(k) are switched. In later revisions of the paper this error has
been corrected.

In both algorithms we must add bounds checking, e.g. the second
algorithm can be modified to:

```C++
// Special leaves which are the product of two primes
for i := pi(y) to pi(sqrt(y))
	if (a / primes[i] < y)
		l := pi(a / primes[i]) + 1;
		limit := pi(min(b / primes[i], y));
		while l <= limit
			// it is a special leaf
			process(primes[i] * primes[l]);
			l := l + 1
```

process(n)
----------

We have found a special leaf, compute it's contribution 
```phi(x / n, i - 1)``` by counting the number of unsieved elements = x / n
after having removed the multiples of the first b primes from the
sieve array. The code below uses the special counters data structure
from Tomás Oliveira e Silva's paper
"Computing pi(x): the combinatorial method", Revista do DETUA, vol. 4,
no. 6, March 2006, pp. 759-768.

```C++
int64_t count = cnt_query(counters, (x / n) - a);
int64_t phi_n = phi[i] + count;

result -= mu[m] * phi_n;
```
