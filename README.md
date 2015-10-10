primecount-backup
=================
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)

The primecount backup version saves intermediate results to a file once
per hour (configurable using ```--backup[=N]```). If your computer crashes
and you restart the same computation primecount will automatically resume
from the backup files.

The backup functionality works quite nicely but the code quality does not
yet meet my standards to be included into the master branch.

### Binaries
Below are the latest precompiled binaries for Windows 64-bit and Linux x86-64.
These binaries are statically linked and require a CPU (2010 or later) which
supports the POPCNT instruction.

* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-backup-2.3-win64.zip">primecount-backup-2.3-win64.zip</a>, 421K
* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-backup-2.3-linux-x64.tar.gz">primecount-backup-2.3-linux-x64.tar.gz</a>, 931K

SHA1 checksums of the files:
```sh
7ee9eeb2bb18bf07b1354deb47dfe0d3486e71c8  primecount-backup-2.3-win64.zip
9543800ac3adf5f2f411842f5774043a251f8176  primecount-backup-2.3-linux-x64.tar.gz
```

### Backup usage example
```sh
# We start a computation and then simulate a crash using Ctrl + C
$ ./primecount 1e22 --S2_hard --status --backup=60

=== S2_hard(x, y) ===
Computation of the hard special leaves
x = 10000000000000000000000
y = 2418290079
z = 4135153217076
c = 6
alpha = 112.247
threads = 8

Status: 43.2%, Load balance: 97%^C
```

```sh
# Now when we rerun the same computation primecount resumes from the backup file
$ ./primecount 1e22 --S2_hard --status --backup=60

=== S2_hard(x, y) ===
Computation of the hard special leaves
x = 10000000000000000000000
y = 2418290079
z = 4135153217076
c = 6
alpha = 112.247
threads = 8

--- Resuming from S2_hard.txt ---
low = 3842113537
segment_size = 2097152
segments_per_thread = 1
S2_hard = 21554940363642178312
Seconds = 3624.303
Status = 33.1%
```

### Command-line options
```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Options:

  -b<N>, --backup=<N>       Store backup file every N minutes (default 60)
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
