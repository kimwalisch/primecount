# primecount MPI

[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primecount.svg)](https://github.com/kimwalisch/primecount/releases)

This is the distributed version of primecount which uses the
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) library
for inter-process communication and which automatically distributes
the computation onto cluster nodes. primecount MPI scales nearly
linearly (up to 96.8% efficiency) for large input values and has been
tested successfully with clusters of up to 50 nodes.

## Build instructions (Unix-like OSes)

First install the prerequisites:
```sh
sudo apt-get install g++ make cmake libopenmpi-dev openmpi-bin
```

Then build primecount MPI using:
```sh
cmake -DWITH_MPI=ON .
make -j
```

## Usage example

```sh
# Distribute pi(10^23) computation onto 30 cluster nodes
mpiexec -n 30 --map-by node ./primecount 1e23 --status
```

```--map-by node``` ensures that only one primecount process will be
created on each cluster node. This is important for performance as
primecount will by default use all available CPU cores on each cluster
node using [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
multi&#8209;threading.

## Performance tips

* You should create only one process per cluster node as primecount MPI
  will by default use all available CPU cores on each cluster node using
  [OpenMP](https://en.wikipedia.org/wiki/OpenMP) multi-threading.
* Be careful when submitting multiple future primecount MPI jobs, they
  might be run simultaneously on the same hardware instead of one after
  the other. This will of course deteriorate performance.
* Ideally all cluster nodes should have an equal number of CPU cores
  because primecount MPI evenly (statically) distributes the work among
  the cluster nodes.
* If possible you should
  [enable transparent huge pages](https://github.com/kimwalisch/primecount#performance-tips)
  on all cluster nodes in order to reduce [TLB (translation lookaside buffer)](https://en.wikipedia.org/wiki/Translation_lookaside_buffer)
  cache misses. This will usually provide a minor speedup.

## Benchmark pi(10<sup>23</sup>)

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
    <td>30,807.46</td>
    <td>1.00 x</td>
    <td>100.00%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>5</td>
    <td>80</td>
    <td>6,645.98</td>
    <td>4.63 x</td>
    <td>92.60%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10</td>
    <td>160</td>
    <td>3,204.08</td>
    <td>9.61 x</td>
    <td>96.10%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>20</td>
    <td>320</td>
    <td>1,590.61</td>
    <td>19.36 x</td>
    <td>96.80%</td>
  </tr>
  </tr>
  <tr align="right">
    <td>30</td>
    <td>480</td>
    <td>1,064.26</td>
    <td>28.94 x</td>
    <td>96.46%</td>
  </tr>
  <tr align="right">
    <td>40</td>
    <td>640</td>
    <td>827.73</td>
    <td>37.21 x</td>
    <td>93.02%</td>
  </tr>
  <tr align="right">
    <td>50</td>
    <td>800</td>
    <td>718.15</td>
    <td>42.89 x</td>
    <td>85.78%</td>
  </tr>
</table>

The pi(10<sup>23</sup>) benchmark above was run on an
[EC2 cluster](https://aws.amazon.com/ec2/) where each cluster node had
2 CPUs of type Intel Xeon E5-2680 v2 (2.80GHz, 8 CPU cores, 16 threads).
The efficiency drops beyond 40 cluster nodes because the input
10<sup>23</sup> is too small for such a large number of nodes.
For 10<sup>24</sup> and 50 cluster nodes the efficiency is 93,65%.

## Command-line options

```
Usage: primecount x [options]
Count the number of primes less than or equal to x (<= 10^31).

Options:

  -d, --deleglise-rivat  Count primes using the Deleglise-Rivat algorithm
  -g, --gourdon          Count primes using Xavier Gourdon's algorithm
  -l, --legendre         Count primes using Legendre's formula
      --lehmer           Count primes using Lehmer's formula
      --lmo              Count primes using Lagarias-Miller-Odlyzko
  -m, --meissel          Count primes using Meissel's formula
      --Li               Approximate pi(x) using the logarithmic integral
      --Li-inverse       Approximate the nth prime using Li^-1(x)
  -n, --nth-prime        Calculate the nth prime
  -p, --primesieve       Count primes using the sieve of Eratosthenes
      --phi <X> <A>      phi(x, a) counts the numbers <= x that are not
                         divisible by any of the first a primes
      --Ri               Approximate pi(x) using Riemann R
      --Ri-inverse       Approximate the nth prime using Ri^-1(x)
  -s, --status[=NUM]     Show computation progress 1%, 2%, 3%, ...
                         Set digits after decimal point: -s1 prints 99.9%
      --test             Run various correctness tests and exit
      --time             Print the time elapsed in seconds
  -t, --threads=NUM      Set the number of threads, 1 <= NUM <= CPU cores.
                         By default primecount uses all available CPU cores.
  -v, --version          Print version and license information
  -h, --help             Print this help menu

Advanced options for the Deleglise-Rivat algorithm:

  -a, --alpha=NUM        Set tuning factor: y = x^(1/3) * alpha
      --P2               Compute the 2nd partial sieve function
      --S1               Compute the ordinary leaves
      --S2-trivial       Compute the trivial special leaves
      --S2-easy          Compute the easy special leaves
      --S2-hard          Compute the hard special leaves
```
