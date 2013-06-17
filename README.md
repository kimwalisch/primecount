primecount
==========

The primecount project contains multiple C++ implementations of the prime counting function. So far I have implemented Legendre's, Meissel's and Lehmer's formulas. All implementations are fully parallelized using OpenMP. primecount can easily be built on any Unix-like operating system (using GNU make) and offers both a command-line program and a C++ library with an intuitive API.

How to build it
---------------
  
First you need to build and install my primesieve library:  
```
$ svn checkout http://primesieve.googlecode.com/svn/trunk/ primesieve
$ cd primesieve
$ make lib
$ sudo make install
$ cd ..
```
  
Then build primecount using:
```
$ make
```
  
How to use it
-------------
  
```
$ bin/./primecount 10^13
```
