Computing special leaves
========================

This file describes how to implement the algorithm to find the
special leaves (sub-algorithm of the Lagarias-Miller-Odlyzko prime
counting algorithm) as described in section 5 (pages 556-557) in the
paper "Computing pi(x) The Meissel-Lehmer Method", Mathematics of
Computation, 44 (1985), pp. 537–560, by J. C. Lagarias, V. S. Miller
and A. M. Odlyzko.

The algorithm described in this document is implemented in
[pi_lmo5.cpp](../src/pi_lmo5.cpp).

Definitions and bounds
----------------------

* pi(x) counts the primes below x
* lpf(n) denotes the smallest prime factor of n
* mu(n) is the Möbius function
* square_free(n) returns true if n is square free
* x^(2/5) > y > x^(1/3)
* j ≤ pi(sqrt(y))
* a = lower bound of the current segment (algorithm uses the segmented sieve of Eratosthenes)
* b = upper bound of the current segment (...)
* k is a small constant e.g. 5, used to pre-sieve the multiples of the first k primes

Prerequisities
--------------

* Generate a pi[n] lookup table of size y
* Generate a lpf[n] lookup table of size y
* Generate a mu[n] lookup table of size y
* Generate a square_free[n] table of size y
* Calculate for each prime pk ≤ sqrt(y) two parallel tables Ak and Mk
* The value of Ak(j) is the jth square-free n ≤ y such that lpf(n) = pk
* The value of Mk(j) is mu(Ak(j))
* For each prime pk ≤ sqrt(y) we store a table called Nk such that l = Nk(j) satisfies Ak(l - 1) < j*pk ≤ Ak(l)

<p>These tables may be computed in time O(y log x) and take up
space O(y log log x).</p>

### Notes
Mk(j) does not seem to be used in the algorithm!?

Ak table
--------

Calculate for each prime pk ≤ sqrt(y) the parallel table Ak. The value
of Ak(j) is the jth square-free n ≤ y such that lpf(n) = pk. Ak(j) can
be calculated efficiently using the code below:

```C++
// A[k][j]
vector<vector<int> > A(/* size = */ pi[sqrt(y)] + 1, vector<int>(1));

for (int n = 1; n <= y; n++)
    if (square_free[n] && lpf[n] <= sqrt(y))
        A[pi[lpf[n]]].push_back(n);
```

### Examples

* ```k = 5```, thus ```primes[5] = 11```
*  11 is the 1st square free n with lpf(n) = 11, thus ```A[5][1] = 11```
* 143 is the 2nd square free n with lpf(n) = 11, thus ```A[5][2] = 143```
* 187 is the 3rd square free n with lpf(n) = 11, thus ```A[5][3] = 187```

Nk table
--------

For each prime pk ≤ sqrt(y) we store a table called Nk such that
l = Nk(j) satisfies Ak(l-1) < j*pk ≤ Ak(l). Our algorithm simply
iterates over all indices of the Ak table and for each Ak(l) finds
a j which satisfies Ak(l-1)/pk < j ≤ Ak(l)/pk. We use a vector of maps
as our data structure, querying N(k, j) uses O(log n) operations.

```C++
// N_maps[k][j]
vector<map<int, int> > N_maps(A.size());

for (int k = 1; k < A.size(); k++) {
    for (int l = 1; l < A[k].size(); l++) {
        int j = A[k][l] / primes[k];
        N_maps[k][j] = l;
    }
}

int N(int k, int j) {
    // find the smallest key >= j using binary search
    auto iter = N_maps[k].lower_bound(j);
    if (iter == N_maps[k].end())
        return INT_MAX;

    // l = N[k][key /* >= j */]
    int l = iter->second;
    return l;
}
```

### Example

* ```k = 5```, ```j = 16```, thus ```j * primes[k] = 176```
* ```l = N(k, j)```, thus ```N(5, 16) = 3```
* ```A[k][l-1] < j * primes[k] <= A[k][l] ```, thus ```A[5][2] < 176 <= A[5][3]```
* ```143 < 176 <= 187```

Algorithm
---------

<p>To find all special n ∈ [a, b) we use these two tables together with
the procedure from section 5 (page 557, figure 2).</p>

```C++
// Special leaves with lpf[n] <= sqrt(y)
for i := k to pi[sqrt(y)]
    for j := i + 1 to pi[sqrt(y)]
        l :=  N(j, (a - 1) / primes[j] + 1)
        while A[j][l] <= b / primes[i]
            // it is a special leaf
            process(primes[i] * A[j][l]);
            l := l + 1

// Special leaves which are the product of two primes
for i := pi[sqrt(y)] to pi[y]
    l := pi[a / primes[i]] + 1;
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

Segmentation bounds (2nd algorithm)
-----------------------------------

In each segment and for each ```primes[i]``` we need to process all
the special leaves that satisfy:

```C++
// special leaf
n = primes[i] * primes[l]
segment_low <= x / n < segment_high
```

The code below shows how to find the lower and upper bound for ```l```:

```C++
for (int i = pi[sqrt(y)]; i + 1 < pi[y]; i++)
{
    int prime = primes[i + 1];
    int l_min = pi[max(x / (prime * segment_high), y / prime)] + 1;
    int l_max = pi[min(x / (prime * segment_low ), y)];

    if (prime >= primes[l_max])
        break;

    for (int l = l_min; l <= l_max; l++)
    {
        n = prime * primes[l];
        process(n);
    }
}
```

### Optimization

The code further up uses slow integer division operations to
calculate the ```l``` bounds. We can speed up the code by
using an ```l_max[]``` array and backwards iterating over ```l```.

```C++
int special_leaf_threshold = max(x / segment_high, y);

for (int i = pi[sqrt(y)]; i + 1 < pi[y]; i++)
{
    int prime = primes[i + 1];
    int l = l_max[i];
    if (prime >= primes[l])
        break;

    for (prime * primes[l] > special_leaf_threshold; l--)
        process(prime * primes[l]);

    l_max[i] = l;
}
```

process(n)
----------

We have found a special leaf, compute it's contribution 
```phi(x / n, i - 1)``` by counting the number of unsieved elements ≤ x / n
after having removed the multiples of the first ```i - 1``` primes from
the sieve array. The code below uses the special counters data structure
from Tomás Oliveira e Silva's paper
"Computing pi(x): the combinatorial method", Revista do DETUA, vol. 4,
no. 6, March 2006, pp. 759-768.

```C++
int64_t count = cnt_query(counters, (x / n) - segment_low);
int64_t phi_xn = phi[i] + count;

result += phi_xn;
```
