primecount MPI
==============
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primecount/blob/master/COPYING)

This is a distributed version of primecount which uses the
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) library for
inter-process communication and which automatically distributes
the computation onto cluster nodes. Breaking the next world record pi(10<sup>28</sup>)
would take about 6 years on the fastest shared memory server currently
available. Hence it has become necessary to go distributed!

Computing pi(10<sup>28</sup>) using primecount MPI requires up to 500 gigabytes
of memory per cluster node! That is a lot, it might take a few years until
such clusters become more widely available. I expect the pi(10<sup>28</sup>)
computation to take about 130 CPU core years using primecount MPI.

Build instructions (Unix-like OSes)
-----------------------------------

First install the prerequisites:
```sh
sudo apt-get install g++ make libopenmpi-dev openmpi-bin
```

Then download
[primecount-3.4.tar.gz](https://dl.bintray.com/kimwalisch/primecount/primecount-3.4.tar.gz)
and build it:
```sh
./build.sh --enable-mpi
```

Usage example
-------------

```sh
# Distribute pi(10^23) computation onto 100 cluster nodes
mpiexec -n 100 -bynode -hostfile my_hosts ./primecount 1e23 --status
```

Note that you should create only one process per cluster node as
primecount MPI will by default use all available CPU cores
on each node using [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
multi-threading.

Benchmark pi(10<sup>23</sup>)
-----------------------------
<table>
  <tr align="center">
    <td><b>Cluster nodes</b></td>
    <td><b>CPU cores</b></td>
    <td><b>Seconds</b></td>
    <td><b>Speedup</b></td>
    <td><b>Efficiency</b></td>
  </tr>
  <tr align="right">
    <td>1</td>
    <td>16</td>
    <td>34,952.48</td>
    <td>1.00 x</td>
    <td>100.00%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>5</td>
    <td>80</td>
    <td>7,968.25</td>
    <td>4.38 x</td>
    <td>87.60%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10</td>
    <td>160</td>
    <td>3,914.31</td>
    <td>8,93 x</td>
    <td>89.30%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>20</td>
    <td>320</td>
    <td>1,922.08</td>
    <td>18.18 x</td>
    <td>90.90%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>30</td>
    <td>480</td>
    <td>1,292.46</td>
    <td>27.04 x</td>
    <td>90.13%</td>
  </tr>
  <tr align="right">
    <td>40</td>
    <td>640</td>
    <td>972.62</td>
    <td>35.93 x</td>
    <td>89.82%</td>
  </tr>
  <tr align="right">
    <td>50</td>
    <td>800</td>
    <td>784.99</td>
    <td>44.52 x</td>
    <td>89.04%</td>
  </tr>
</table>

The pi(10<sup>23</sup>) benchmark above was run on an
[EC2 cluster](https://aws.amazon.com/ec2/) where each cluster node had
2 CPUs of type Intel Xeon E5-2680 v2 (2.80GHz, 8 CPU cores, 16 threads).
The cluster has been set up using the
[StarCluster](http://star.mit.edu/cluster/) tool. The efficiency drops
slightly beyond 30 cluster nodes, the author thinks this is because
the input 10<sup>23</sup> is too small.

Command-line options
--------------------

```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Options:

  -d,    --deleglise_rivat  Count primes using Deleglise-Rivat algorithm
         --legendre         Count primes using Legendre's formula
         --lehmer           Count primes using Lehmer's formula
  -l,    --lmo              Count primes using Lagarias-Miller-Odlyzko
  -m,    --meissel          Count primes using Meissel's formula
         --Li               Approximate pi(x) using the logarithmic integral
         --Li_inverse       Approximate the nth prime using Li^-1(x)
  -n,    --nthprime         Calculate the nth prime
  -p,    --primesieve       Count primes using the sieve of Eratosthenes
  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...
                            [N] digits after decimal point e.g. N=1, 99.9%
         --test             Run various correctness tests and exit
         --time             Print the time elapsed in seconds
  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores
  -v,    --version          Print version and license information
  -h,    --help             Print this help menu

Advanced Deleglise-Rivat options:

  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)
         --P2               Only compute the 2nd partial sieve function
         --S1               Only compute the ordinary leaves
         --S2_trivial       Only compute the trivial special leaves
         --S2_easy          Only compute the easy special leaves
         --S2_hard          Only compute the hard special leaves
```
