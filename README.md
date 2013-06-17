primecount
==========

The primecount project contains multiple C++ implementations of the prime counting function. So far I have implemented Legendre's, Meissel's and Lehmer's formulas. All implementations are fully parallelized using OpenMP. primecount can easily be built on any Unix-like operating system (using GNU make) and offers both a command-line program and a C++ library with an intuitive API.

### How to build it
primecount depends on the author's primesieve libary (version 4.3 or later). To download, build and install the latest primesieve and libprimesieve version run:
```
$ sh install_primesieve.sh
```
Then build primecount using GNU make and optionally install it:
```
$ make
$ sudo make install
```

### How to use it
```
$ bin/./primecount 10^13
```
