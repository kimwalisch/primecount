# primecount-backup

[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primecount.svg)](https://github.com/kimwalisch/primecount/releases)

The primecount backup version saves intermediate results to a backup file.
If your computer crashes or if you interrupt a computation you can resume
the same computation from the backup file. For pi(x) computations that
take weeks or even months to compute the backup functionality is very
important. David Baugh and myself have used primecount-backup to compute
pi(10^27) and many other records.

The primecount backup version generates the following files:

<table>
    <tr>
        <td><code>primecount.backup</code></td>
        <td>Backup file with intermediate results</td>
    </tr>
    <tr>
        <td><code>primecount.log</code></td>
        <td>Detailed log file</td>
    </tr>
    <tr>
        <td><code>results.txt</code></td>
        <td>Log file with final results only</td>
    </tr>
</table>

## Build instructions

You need to have installed a C++ compiler and CMake. Ideally
primecount should be compiled using a C++ compiler that supports both
OpenMP and 128-bit integers (e.g. GCC, Clang, Intel C++ Compiler).

```sh
cmake .
make -j
```

* [Detailed build instructions](doc/BUILD.md)

## Binaries

Below are the latest precompiled primecount binaries for
Windows and Linux. These binaries are statically linked and
require a CPU which supports the POPCNT instruction (2008 or
later).

* [primecount-backup-5.2-win64.zip](https://github.com/kimwalisch/primecount/releases/download/v5.2-backup/primecount-backup-5.2-win64.zip), 626 KB
* [primecount-backup-5.2-linux-x64.tar.xz](https://github.com/kimwalisch/primecount/releases/download/v5.2-backup/primecount-backup-5.2-linux-x64.tar.xz), 916 KB
* [primecount-backup-5.2-macOS-x64.zip](https://github.com/kimwalisch/primecount/releases/download/v5.2-backup/primecount-backup-5.2-macOS-x64.zip), 427 KB

## Backup usage example

```sh
# We start a computation and interrupt it using Ctrl + C
$ ./primecount 1e21 --B -s

=== B(x, y) ===
x = 1000000000000000000000
y = 351720000
alpha_y = 35.172
threads = 8

Status: 38.7%^C
```

```sh
# Now we resume the computation from the backup file
$ ./primecount --resume

=== B(x, y) ===
x = 1000000000000000000000
y = 351720000
alpha_y = 35.172
threads = 8

Resuming from primecount.backup
MD5 checksum: OK
Status: 38.7%
```

The ```primecount.backup``` file is updated using atomic writes in order to
prevent crashes from corrupting the ```primecount.backup``` file. We also
calculate an MD5 checksum of the  ```primecount.backup``` content before writing
to disk and verify that checksum when resuming in order to protect from
harddisk errors.

Note that the ```primecount.backup``` file is not tied to the PC on which the
computation was originally started. You can safely copy the
```primecount.backup``` file to another PC and resume the computation there.
If the new PC has a different number of CPU cores primecount will by default
resume the computation using all available CPU cores (unless you have
specified the number of threads using ```--threads=NUM```).

## Batch processing

It is possible to create a ```worktodo.txt``` file with a list of
numbers to compute e.g.:

```sh
# worktodo.txt
10000000
1e15
1e15 --alpha-y=10 --threads=4
1e14 --AC
1e18 --D
```

Then you can process all numbers from ```worktodo.txt``` using:

```sh
$ scripts/worktodo.sh
```

The results will be stored in ```results.txt``` and extended
details are logged into ```primecount.log```.

## Verifying pi(x) results

Record pi(x) computations may take many months to complete. For such long
running computations it is important to double check the pi(x) computation in
order to protect from hardware and software errors. The first thing you can do
to reduce the risk of a pi(x) miscalculation due to hardware errors is using
[ECC memory](https://en.wikipedia.org/wiki/ECC_memory).

In order to double check and verify a pi(x) computation you have to run the
same pi(x) computation a second time but this time you manually specify a slightly
different ```alpha_y``` or ```alpha_z``` tuning factor (using e.g. ```--alpha-y=NUM```).
Doing this the results of the many formulas of Gourdon's algorithm will be
completely different from the first run but if the pi(x) results of the 1st and
2nd run match then the computation has been verified successfully!

## Command-line options

```
Usage: primecount x [options]
Count the number of primes less than or equal to x (<= 10^31).

Backup options:

  -b, --backup=FILENAME    Set the backup filename. The default backup
                           filename is primecount.backup.

  -r, --resume[=FILENAME]  Resume the last computation from the
                           primecount.backup file. If another backup
                           filename is provided the computation is resumed
                           from that backup file.
Options:

  -g, --gourdon            Count primes using Xavier Gourdon's algorithm.
                           This is the default algorithm.
  -l, --legendre           Count primes using Legendre's formula
  -m, --meissel            Count primes using Meissel's formula
      --Li                 Approximate pi(x) using the logarithmic integral
      --Li-inverse         Approximate the nth prime using Li^-1(x)
  -n, --nth-prime          Calculate the nth prime
  -p, --primesieve         Count primes using the sieve of Eratosthenes
      --phi <X> <A>        phi(x, a) counts the numbers <= x that are not
                           divisible by any of the first a primes
      --Ri                 Approximate pi(x) using Riemann R
      --Ri-inverse         Approximate the nth prime using Ri^-1(x)
  -s, --status[=NUM]       Show computation progress 1%, 2%, 3%, ...
                           Set digits after decimal point: -s1 prints 99.9%
      --test               Run various correctness tests and exit
      --time               Print the time elapsed in seconds
  -t, --threads=NUM        Set the number of threads, 1 <= NUM <= CPU cores.
                           By default primecount uses all available CPU cores.
  -v, --version            Print version and license information
  -h, --help               Print this help menu

Advanced options for Xavier Gourdon's algorithm:

      --alpha-y=NUM        Set tuning factor: y = x^(1/3) * alpha_y
      --alpha-z=NUM        Set tuning factor: z = y * alpha_z
      --AC                 Compute the A + C formulas
      --B                  Compute the B formula
      --D                  Compute the D formula
      --Phi0               Compute the Phi0 formula
      --Sigma              Compute the 7 Sigma formulas
```
