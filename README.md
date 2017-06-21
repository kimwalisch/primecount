# primecount-backup

[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primecount/blob/master/COPYING)

The primecount backup version saves intermediate results to a backup file.
If your computer crashes or if you interrupt a computation you can resume
the same computation from the backup file. For pi(x) computations that
take weeks or even months to compute the backup functionality is very
important. David Baugh and myself have used primecount-backup to compute
pi(10^27) and many other records.

## Binaries

Below are the latest precompiled primecount binaries for
Windows and Linux. These binaries are statically linked and
require a CPU which supports the POPCNT instruction (2010 or
later).

* [primecount-3.7-win64.zip](https://github.com/kimwalisch/primecount/releases/download/v3.7/primecount-backup-3.7-win64.zip), 479 KB
* [primecount-3.7-linux-x64.tar.gz](https://github.com/kimwalisch/primecount/releases/download/v3.7/primecount-backup-3.7-linux-x64.tar.gz), 1.07 MB

## Backup usage example

```sh
# We start a computation and interrupt it using Ctrl + C
$ ./primecount 1e19 --P2 -s

=== P2(x, y) ===
Computation of the 2nd partial sieve function
x = 10000000000000000000
y = 75434780
z = 132564846082
alpha = 35.014
threads = 8

Status: 52%^C
```

```sh
# Now we resume the computation from the backup file
$ ./primecount --resume

=== P2(x, y) ===
Computation of the 2nd partial sieve function
x = 10000000000000000000
y = 75434780
z = 132564846082
alpha = 35.014
threads = 8

=== Resuming from primecount.backup ===
low = 68652367874
thread_distance = 7989059776
pix_total = 2871709096
p2 = 56268889109074621
Seconds: 10.153

Status: 52%
```

## Command-line options

```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Backup options:

  -b, --backup=<filename>   Set the backup filename. The default backup
                            filename is primecount.backup

  -r, --resume[=filename]   Resume the last computation from the
                            primecount.backup file. If another backup
                            filename is provided the computation is resumed
                            from that backup file
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
         --phi=<a>          phi(x, a) counts the numbers <= x that are
                            not divisible by any of the first a primes
         --Ri               Approximate pi(x) using Riemann R
         --Ri_inverse       Approximate the nth prime using Ri^-1(x)
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

Examples:

  primecount 1e13
  primecount 1e13 --nthprime --threads=4
```

## Batch processing

It is possible to create a ```worktodo.txt``` file with a list of
numbers to compute e.g.:

```sh
# worktodo.txt
10000000
1e15
1e15 --alpha=10 --threads=4
1e14 --P2
1e18 --S2_hard
```

Then you can batch process all the numbers from ```worktodo.txt``` using:

```sh
$ sh scripts/worktodo.sh
```

The results will be stored in ```results.txt``` and more details are
logged into ```primecount.log```.
