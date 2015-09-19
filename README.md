primecount with backup functionality
====================================
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)

This primecount version saves intermediate results of P2, S2_easy and S2_hard
to a file once per hour. If your computer crashes and you restart the same
computation primecount will automatically resume from the backup files.

### Binaries
Below are the latest precompiled binaries for Windows 64-bit and Linux x86-64.
These binaries are statically linked and require a CPU (2010 or later) which
supports the POPCNT instruction.

* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-backup-2.2-win64.zip">primecount-backup-2.2-win64.zip</a>, 380K
* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-backup-2.2-linux-x64.tar.gz">primecount-backup-2.2-linux-x64.tar.gz</a>, 936K

SHA1 checksums of the files:
```sh
1815630ea53500e1a00fc42bb940ccf51e8bf380  primecount-backup-2.2-win64.zip
6320ed9cea03f31b1d923340283d6fe3a53fe048  primecount-backup-2.2-linux-x64.tar.gz
```

### Build instructions (Unix-like OSes)
To build primecount-backup you need to have installed a C++ compiler,
GNU make and the GNU Build System (a.k.a. Autotools). To install the
GNU Build System install
[GNU&#160;Autoconf](http://www.gnu.org/software/autoconf/),
[GNU&#160;Automake](http://www.gnu.org/software/automake/) and
[GNU&#160;Libtool](http://www.gnu.org/software/libtool/)
using your package manager.

primecount depends on the author's primesieve library, download it from
http://primesieve.org/downloads and install it using:
```sh
$ ./configure
$ make
$ sudo make install
```

If you are not using Linux then you need to export these variables:
```sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
```

Finally download the latest
[primecount-backup.zip](https://github.com/kimwalisch/primecount/archive/backup.zip)
archive from GitHub and build it using:
```sh
$ ./autogen.sh
$ ./configure
$ make
$ sudo make install
```

If your CPU supports the
[POPCNT instruction](http://en.wikipedia.org/wiki/SSE4#POPCNT_and_LZCNT)
then it is enabled in the build process. POPCNT speeds up primecount by
about 10 percent. If you need maximum portability you can disable POPCNT:
```sh
$ ./configure --disable-popcnt
```

### Build instructions (Microsoft Visual C++)
To build primecount simply open a Visual Studio Command Prompt and execute:
```sh
> nmake -f Makefile.msvc
```

### Command-line options
```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using the prime counting function,
by default the Deleglise-Rivat algorithm (-d) is used.

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
  -s,    --status           Print status info during computation
         --test             Run various correctness tests and exit
         --time             Print the time elapsed in seconds
  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores
  -v,    --version          Print version and license information
  -h,    --help             Print this help menu

Advanced Deleglise-Rivat options:

  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)
         --p2               Only compute the 2nd partial sieve function
         --s1               Only compute the ordinary leaves
         --s2_trivial       Only compute the trivial special leaves
         --s2_easy          Only compute the easy special leaves
         --s2_hard          Only compute the hard special leaves
```
