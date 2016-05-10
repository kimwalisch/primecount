Hacking on primecount
=====================

### Hacking on Unix-like OSes

Clone or fork primecount:
```sh
$ git clone git://github.com/kimwalisch/primecount.git
```

In order to hack primecount you need to have installed a C++ compiler
and the GNU Build System (a.k.a. Autotools). To install the GNU Build
System install
[GNU Autoconf](http://www.gnu.org/software/autoconf/),
[GNU Automake](http://www.gnu.org/software/automake/) and
[GNU Libtool](http://www.gnu.org/software/libtool/) using your packet
manager.

Build primecount using:
```sh
$ ./build.sh
```

### Adding a new prime counting function implementation

* Add new function signature e.g. ```pi_lmo(x)``` to [include/primecount.hpp](include/primecount-internal.hpp)
* Create source file ```src/lmo/pi_lmo.cpp``` with the new function
* Add ```src/lmo/pi_lmo.cpp``` to [Makefile.am](Makefile.am)
* Add ```src\lmo\pi_lmo.obj``` to [Makefile.msvc](Makefile.msvc)
* Add ```pi_lmo(x)``` to [src/app/main.cpp](src/app/main.cpp)
* Add ```pi_lmo(x)``` to ```bool test()``` in [src/test.cpp](src/test.cpp)
* Add ```--lmo``` command-line option summary to [src/app/help.cpp](src/app/help.cpp)
* Add ```OPTION_LMO``` to [src/app/cmdoptions.hpp](src/app/cmdoptions.hpp)
* Add ```OPTION_LMO``` to optionMap in [src/app/cmdoptions.cpp](src/app/cmdoptions.cpp)

### Versioning

* Increase version number in [include/primecount.hpp](include/primecount.hpp)
* Increase version number in _**Build instructions**_ section in [README.md](README.md)
* Increase version number in _**Build instructions**_ section in [doc/primecount-MPI.md](doc/primecount-MPI.md)
* Increase version number in [configure.ac](configure.ac) in ```AC_INIT```
* [Increase Libtool version](http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html) number in [configure.ac](configure.ac) in ```AC_SUBST```
* Update to current year in [src/app/help.cpp](src/help.cpp)

### Release process

* Run tests using ```make check```
* Increase version number (see <a href="#versioning">Versioning</a>)
* Build statically linked primecout binaries and upload them to [https://bintray.com/kimwalisch/primecount](https://bintray.com/kimwalisch/primecount)
* Update _**Precompiled binaries**_ section in [README.md](README.md)
* Update [ChangeLog](ChangeLog)
* Tag the new release in git
* Create new release tarball using ```make dist``` and upload it to [https://bintray.com/kimwalisch/primecount](https://bintray.com/kimwalisch/primecount)
