# Computation of the easy special leaves

In the combinatorial prime counting algorithms the computation of the special leaves is the computationally most expensive task.
In order to speed up that computation Lagarias-Miller-Odlyzko [[1]](#references) have split up the special leaves into easy
special leaves and hard special leaves. The contribution of each easy special leaf can be computed in $O(1)$ using a
```PrimePi[n]``` lookup table whereas the contribution of each hard special leaf requires evaluating the partial sieve function
$\phi(x, a)$ and cannot be computed in $O(1)$.

In the Deleglise-Rivat [[2]](#references) and Gourdon [[3]](#references) prime counting algorithms (which are based on the
Lagarias-Miller-Odlyzko algorithm) the computation of the easy special leaves requires looking up the number of primes below n
with $n < \sqrt{x}$. Since a ```PrimePi[n]``` lookup table of size $\sqrt{x}$ is much too large to be practical, Deleglise-Rivat
[[2]](#references) have suggested segmenting the Interval $[0, \sqrt{x}$[ using a segment size of y (~ $\sqrt[3]{x}\ \log^{3}{x}$).
So instead of using a ```PrimePi[n]``` lookup table of size $\sqrt{x}$ we now use a ```SegmentedPrimePi[n]``` lookup table of size
y which also returns the number of primes ≤ n but requires n to be within the current segment [low, low + y[. This approach was
used in primecount up to version 6.4. However this segment size causes severe scaling issues for large computations > $10^{22}$ as the
```SegmentedPrimePi[n]``` lookup table becomes exceedingly large e.g. at $10^{30}$ its size was 137 GiB in primecount. For this reason
Xavier Gourdon [[3]](#references) suggested using a smaller segment size of $\sqrt{\frac{x}{y}}$ which is orders of magnitude
smaller and generally a good practical improvement.

Here are links to primecount's [PiTable](https://github.com/kimwalisch/primecount/blob/master/src/PiTable.cpp) and
[SegmentedPiTable](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/SegmentedPiTable.cpp) implementations.

# Improving the cache efficiency

The ```SegmentedPrimePi[n]``` lookup table is accessed very frequently in the computation of the easy special leaves (about once for each
easy special leaf) and these memory accesses are non sequential. It is therefore important that the ```SegmentedPrimePi[n]``` fits into
the CPU's fast cache memory. While Xavier Gourdon's smaller segment size is already considerably smaller it is still too large for new
record computations. For this reason I suggest using an **even smaller segment size of $\sqrt[4]{x}$** for the computation of the easy
special leaves. With a segment size of $\sqrt[4]{x}$ the ```SegmentedPrimePi[n]``` lookup table fits into the CPU's cache even for record
computations e.g. at $10^{30}$ the ```SegmentedPrimePi[n]``` is only about 2 MiB in primecount. A segment
size of $\sqrt[4]{x}$ does not deteriorate the runtime complexity of the algorithm because the segmented sieve of Eratosthenes which is
used to initialize the ```SegmentedPrimePi[n]``` lookup table has the same runtime complexity as the sieve of Eratosthenes as long as
the segment size is not smaller than the square root of the total sieving distance (which is $\sqrt{x}$).

Note that Deleglise-Rivat [[2]](#references) have split up the easy special leaves into many formulas and suggest using segmentation only for the 2
formulas that need to lookup the number of primes < $\sqrt{x}$, whereas all other formulas that only need to lookup the number of
primes ≤ y should be computed without segmentation. As a ```PrimePi[n]``` lookup table of size y is much too large to fit into the CPU's
cache and as the ```PrimePi[n]``` lookup table is accessed in random order, I suggest segmenting all easy special leaves formulas that
are computationally expensive using a segment size of $\sqrt[4]{x}$ in order to improve performance. To reduce the amount of work for
the programmer it is best to sieve the interval $[0, \sqrt{x}[$ only once and compute all easy special leaf formulas within that sieve.

Extra care needs to be used when segmenting the formulas that compute consecutive identical easy leaves more efficiently, sometimes these
leaves are named clustered easy leaves [[4]](#references). In the Deleglise-Rivat algorithm the $W_3$ and $W_5$ formulas compute clustered easy
leaves. These formulas need to access ```PrimePi[n]``` values with n ≤ y, but some of these memory accesses (i.e. those that compute how
many consecutive leaves are identical) may be outside of the segment [low, low + segment_size[. For these memory accesses I suggest using
a ```PrimePi[n]``` lookup table of size y instead of the ```SegmentedPrimePi[n]``` lookup table. Note that it is important for performance
to segment the clustered easy leaves as there is a proportionally large number of these leaves and their computation is expensive.

 # Parallel computation and load-balancing

So far we have focused on improving the cache efficiency of the computation of the easy special leaves. Now we will have a look at
how to parallelize the computation of the easy special leaves so that the algorithm scales well. Generally parallel algorithms
scale well on current CPU architectures if they accomplish the 3 properties below:

* Each thread only operates on its own tiny chunk of memory that fits into the CPU's cache.
* All threads must be independent from each other (i.e. require no synchronization).
* The work must be distributed evenly among all threads in order to avoid load imbalance.

A segment size of $\sqrt[4]{x}$ already accomplishes the first property. So next we have to design our parallel algorithm in a way that
all threads are independent from each other. Luckily Xavier Gourdon [[3]](#references) already devised an idea for how to do this: **at the start of
each new segment [low, low + segment_size[ each thread computes ```PrimePi(low)``` using a prime counting function implementation**
in $O(low^{\frac{2}{3}})$ or less. The result of ```PrimePi(low)``` is required to initialize the ```SegmentedPrimePi[n]``` lookup table
for the current segment [low, low + segment_size[. This algorithm has been implemented in primecount-7.0
(see [SegmentedPiTable.cpp](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/SegmentedPiTable.cpp), [AC.cpp](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/AC.cpp)), it improved performance
by more than 2x at $10^{23}$ on my dual-socket AMD EPYC server compared to primecount-6.4 which used a larger segment size and
required frequent synchronization of threads. It is important to ensure that the additional pre-computations do not deteriorate
the runtime complexity of the algorithm. When sieving up to $\sqrt{x}$ using a segment size of $\sqrt[4]{x}$ there will by exactly $\sqrt[4]{x}$
segments. For each segment we need to compute ```PrimePi(low)``` with low < $\sqrt{x}$. Hence in total the additional pre-computations
have a runtime complexity of $O(\sqrt{x}^{\frac{2}{3}}\sqrt[4]{x}) = O(x^{\frac{7}{12}})$ which does not deteriorate the overall runtime complexity
of the algorithm.

Lastly we have to ensure that the work is distributed evenly among all threads. The easy special leaves are distributed very
unevenly, most of the leaves are located below y (~ $\sqrt[3]{x}\ \log^{3}{x}$) whereas above y the number of leaves slowly decreases and
they become more and more sparse as they approach $\sqrt{x}$. Hence it is critical that the region below y is distributed evenly
among all threads. Based on my benchmarks **a small segment size of $\sqrt[4]{x}$ evenly distributes the work** even on servers with a
large number of CPU cores such as my dual-socket AMD EPYC server with 196 threads. Using a segment size larger than $\sqrt[4]{x}$ such
as $\sqrt[3]{x}$ or y causes significant load imbalance (i.e. some threads will be assigned much more work than others and keep on
computing after most of the threads have already finished their computations) which severely deteriorates performance especially
on PCs and servers with a large number of CPU cores. Above y there are much fewer easy special leaves hence the segment size can
be increased by a small constant factor (16 in primecount) in order to reduce the pre-computation overhead, provided that the
new segment size still fits into the CPU's cache.

# References

1. J. C. Lagarias, V. S. Miller, and A. M. Odlyzko, Computing pi(x): The Meissel-Lehmer method, Mathematics of Computation, 44 (1985), pp. 537–560.
2. M. Deleglise and J. Rivat, "Computing pi(x): The Meissel, Lehmer, Lagarias, Miller, Odlyzko Method", Mathematics of Computation, Volume 65, Number 213, 1996, pp 235–245.
3. Xavier Gourdon, Computation of pi(x) : improvements to the Meissel, Lehmer, Lagarias, Miller, Odllyzko, Deléglise and Rivat method, February 15, 2001.
4. Tomás Oliveira e Silva, Computing pi(x): the combinatorial method, Revista do DETUA, vol. 4, no. 6, March 2006, pp. 759-768.
