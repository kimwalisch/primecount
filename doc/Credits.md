# Credits

My thanks to the following people, without whom primecount wouldn't be what it is today. (Please let me know if I've mistakenly omitted anyone.)

* **David Baugh**<br/>
David was the most extensive user of primecount during early
development. He reported numerous bugs and suggested a lot
of useful features e.g. backup functionality & standalone
flags for all deleglise-rivat formulas. David also broke many
world records using primecount e.g. [A006880](https://oeis.org/A006880),
[A006988](https://oeis.org/A006988), [A122121](https://oeis.org/A122121),
[A040014](https://oeis.org/A040014).

* [Dana Jacobsen](https://github.com/danaj)<br/>
Dana explained to me in great detail the optimizations used in his
phi(x, a) implementation and also the implementation of the 
inverse logarithmic integral using Newton's method.

* [Curtis Seizert](https://github.com/curtisseizert)<br/>
Curtis partially implemented Xavier Gourdon's combinatorial
prime counting algorithm and wrote
[a paper](https://github.com/curtisseizert/CUDApix/blob/master/Deconvoluting%20Deleglise-Rivat.pdf)
comparing the Deleglise-Rivat algorithm, Tomás Oliveira e Silva's
version of the Deleglise-Rivat algorithm and Gourdon's algorithm.
This improved my understanding of Gourdon's algorithm.

* **Douglas Staple**<br/>
After computing pi(10^26) using his own program Douglas Staple
published [a paper](https://arxiv.org/pdf/1503.01839.pdf) with
his improvements to the Deleglise-Rivat algorithm. I use Douglas
Staple's improvements for the computation of the ordinary
leaves and the multi-threading of the hard special leaves.

* **Christian Bau**<br/>
Christian Bau wrote the first open source implementation of the
extended Meissel-Lehmer prime counting algorithm in 2003. I use
Christian's FactorTable idea to reduce the memory usage and
also his fast integer division trick (use 32-bit instead of 64-bit
whenever possible).

* [Dennis Mitchell](https://codegolf.stackexchange.com/users/12012/dennis)<br/>
Dennis wrote a highly optimized
[implementation](https://codegolf.stackexchange.com/a/74372/52196)
(in 2016) of the meissel prime counting algorithm which was
faster than primecount for x <= 10<sup>10</sup>. Dennis's
implementation replaced integer division by multiplication & bit
shifts. After studying his implementation I also implemented this
trick in primecount.

* [ridiculousfish](https://github.com/ridiculousfish)<br/>
primecount uses ridiculousfish's
[libdivide library](https://github.com/ridiculousfish/libdivide)
to replace integer division by multiplication & bit shifts. This
optimization speeds up the computation of the easy special leaves
by about 40%. ridiculousfish added a branchfree divider to libdivide
specifically for primecount.
