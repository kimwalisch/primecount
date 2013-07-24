primecount
==========
primecount is a command-line program and C++ library that counts the number of primes below an integer x < 2^63. It does so by using an efficient implementation of the prime counting function pi(x) which is orders of magnitude faster than counting primes using the sieve of Eratosthenes. So far primecount offers the option to count primes using Legendre's, Meissel's and Lehmer's formulas. All pi(x) implementations are fully parallelized using OpenMP.

### Formulas

{|
|  Legendre's Formula  ||  ![Legendre's Formula](images/pi_legendre.png)
| -
|  Meissel's Formula   ||  ![Meissel's Formula](images/pi_meissel.png)
| -
|  Lehmer's Formula    ||  ![Lehmer's Formula](images/pi_lehmer.png)
|}

### How to build it
primecount depends on the author's primesieve libary (version 4.3 or later). To download, build and install the latest primesieve version on a Unix-like operating system run:
```
$ sh install_primesieve.sh
```
To build and install primecount using GNU make and the default `c++' compiler run:
```
$ make
$ sudo make install
```

### Usage Examples
```
$ primecount 10^13
$ primecount 10^14 --meissel --threads=2
$ primecount 78498 --nthprime
```
