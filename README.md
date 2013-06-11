primecount
==========

C++ implementations of various prime counting functions (so far Legendre and Meissel).  
This project has been started in June 2013 and is work in progress.

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
