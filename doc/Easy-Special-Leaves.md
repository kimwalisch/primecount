# Computation of the easy special leaves

In the combinatorial prime counting algorithms the computation of the special leaves is the computationally most expensive task.
In order to speed up that computation Deleglise and Rivat have split up the special leaves into easy special leaves and hard
special leaves. The contribution of each easy special leaf can be computed in O(1) using a ```PrimePi[n]``` lookup table whereas the
contribution of each hard special leaf requires evaluating the partial sieve function phi(x, a) and cannot be computed in O(1).

In the computation of the easy special leaves we need to look up the number of primes below n with n < x^(1/2). Since a
```PrimePi[n]``` lookup table of size x^(1/2) is much too large to be practical, Deleglise-Rivat have suggested segmenting the Interval
[0, x^(1/2)[ using a segment size of y (~ x^(1/3) * log(x)^3). So instead of using a ```PrimePi[n]``` lookup of size x^(1/2) we now use
a ```SegmentedPrimePi[n]``` lookup table of size y which also returns the number of primes ≤ n but requires n to be within the current segment
[low, low + y[. This approach was used in primecount up to version 6.4. However this segment size causes severe scaling issues
for large computations > 10^22 as the ```SegmentedPrimePi[n]``` becomes exceedingly large e.g. at 10^30 its size was 137 GiB in
primecount. For this reason Xavier Gourdon suggested using a smaller segment size of sqrt(x/y) which is orders of magnitude
smaller and generally a good practical improvement.

Hare are links to primecount's [PiTable](https://github.com/kimwalisch/primecount/blob/master/src/PiTable.cpp) and
[SegmentedPiTable](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/SegmentedPiTable.cpp) implementations.

# Improving the cache efficiency

The ```SegmentedPrimePi[n]``` lookup table is accessed very frequently in the computation of the easy special leaves (about once for each
easy special leaf) and these memory accesses are non sequential. It is therefore important that the ```SegmentedPrimePi[n]``` fits into
the CPU's cache. While Xavier Gourdon's smaller segment size is already considerably smaller it is still too large for new record
computations. For this reason I suggest using an **even smaller segment size of x^(1/4)** for the computation of the easy special
leaves. With a segment size of x^(1/4) the ```SegmentedPrimePi[n]``` lookup table fits
into the CPU's cache even for record computations e.g. at 10^30 the ```SegmentedPrimePi[n]``` is only about 2 MiB in primecount. A segment
size of x^(1/4) does not deteriorate the runtime complexity of the algorithm because the segmented sieve of Eratosthenes which is
used to initialize the ```SegmentedPrimePi[n]``` lookup table has the same runtime complexity as the sieve of Eratosthenes as long as
the segment size is not smaller than the square root of the total sieving distance.

Note that Deleglise-Rivat have split up the easy special leaves into many formulas and suggest using segmentation only for the 2
formulas that need to lookup the number of primes < x^(1/2), whereas all other formulas that only need to lookup the number of
primes ≤ y should be computed without segmentation. As a ```PrimePi[n]``` lookup table of size y is much too large to fit into the CPU's
cache and as the ```PrimePi[n]``` lookup table is accessed in random order, I suggest segmenting all easy special leaves formulas that
are computationally expensive using a segment size of x^(1/4) in order to improve performance. However special care needs to be
used for the formulas that compute identical consecutive easy leaves more efficiently, sometimes these formulas are named clustered
easy leaves. In the Deleglise-Rivat algorithm the W3 and W5 formulas compute clustered easy leaves. These formulas
need to access ```PrimePi[n]``` values with n ≤ y but some of these memory accesses (i.e. those that compute how many consecutive leaves
are identical) may be outside of the segment [low, low + segment_size[. For these memory accesses I suggest using a ```PrimePi[n]``` lookup
table of size y instead of the ```SegmentedPrimePi[n]``` lookup table.

 # Parallel computation and load-balancing

So far we have focused on improving the cache efficiency of the computation of the easy special leaves. Now we will have a look at
how to parallelize the computation of the easy special leaves so that it scales well. Generally parallel algorithms scale well on
current CPU architectures if they accomplish the 3 properties below:

* Each thread only operates on his own tiny chunk of memory that fits into the CPU's cache.
* All threads must be independent from each other (i.e. require no synchronization).
* The work must be equally shared between all threads in order to avoid load imbalance. 

A segment size of x^(1/4) already accomplishes the first property. So next we have to design our parallel algorithm in a way that
all threads are independent from each other. Luckily Xavier Gourdon already devised an idea for how to do this: **at the start of
each new segment [low, low + segment_size[ each thread computes ```PrimePi[low]``` using a prime counting function implementation in
O(low^(2/3)) or less**. This way no thread requires any data from another thread. This algorithm has been implemented in
primecount-6.5 (see [AC.cpp](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/AC.cpp) &
[SegmentedPiTable.cpp](https://github.com/kimwalisch/primecount/blob/master/src/gourdon/SegmentedPiTable.cpp)), it improved performance
by more than 2x at 10^23 on my dual-socket AMD EPYC server compared to primecount-6.4 which used a larger segment size and
required frequent synchronization of threads. It is important to ensure that the additional pre-computations do not deteriorate
the runtime complexity of the algorithm. When sieving up to x^(1/2) using a segment size of x^(1/4) there will by exactly x^(1/4)
segments. For each segment we need to compute ```PrimePi[low]``` with low < x^(1/2). Hence in total the additional pre-computations
have a runtime complexity of O((x^(1/2))^(2/3) * x^(1/4)) = O(x^(7/12)) which does not deteriorate the overall runtime complexity
of the algorithm.

Lastly we have to ensure that the work is distributed evenly amongst all threads. Most of the easy special leaves are below y
(~ x^(1/3) * log(x)^3), hence it is critical that this region is distributed evenly amongst all threads. Based on my benchmarks
**a segment size of x^(1/4) evenly distributes the work** even on servers with a large number of CPU cores such as my
dual-socket AMD EPYC server with 196 threads. Using a segment size larger than x^(1/4) causes significant load imbalance
which deteriorates performance i.e. some threads will be assigned much more work than others and keep on computing while most
threads have already finished their computations. Above y there are much fewer easy special leaves hence the segment size can be
increased by a small constant factor (16 in primecount) in order to reduce the pre-computation overhead, provided that the new
segment size still fits into the CPU's cache.
