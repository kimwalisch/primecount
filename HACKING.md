Hacking on primecount
=====================

### Hacking on Unix-like OSes

Clone or fork primecount:
```sh
$ git clone https://github.com/kimwalisch/primecount.git
```

In order to hack primecount on Unix-like OSes you need to have
installed a C++ compiler and the GNU Build System (a.k.a Autotools).
To install the GNU Build System you need to install
[GNU Autoconf](http://www.gnu.org/software/autoconf/),
[GNU Automake](http://www.gnu.org/software/automake/) and
[GNU Libtool](http://www.gnu.org/software/libtool/) using your packet
manager.

Generate configure script (only once):
```sh
$ ./autogen.sh
```

Then build primecount using:
```sh
$ ./configure
$ make
```

### Adding a new prime counting function implementation

* Add new function signature e.g. ```pi_lmo(x)``` to [include/primecount.hpp](include/primecount.hpp)
* Create source file ```src/pi_lmo.cpp``` with the new function
* Add ```src/pi_lmo.cpp``` to [Makefile.am](Makefile.am)
* Add ```src\pi_lmo.cpp``` to [Makefile.msvc](Makefile.msvc)
* Add ```OPTION_LMO``` to [src/cmdoptions.hpp](src/cmdoptions.hpp)
* Add ```OPTION_LMO``` to optionMap in [src/cmdoptions.cpp](src/cmdoptions.cpp)
* Add ```pi_lmo(x)``` to [src/primecount.cpp](src/primecount.cpp)
* Add ```pi_lmo(x)``` to ```bool test()``` in [src/test.cpp](src/test.cpp)
* Add ```--lmo``` command-line option summary to [src/help.cpp](src/help.cpp)

### Versioning

* Increase version number in [include/primecount.hpp](include/primecount.hpp)
* Increase version number in _**Build instructions**_ section in [README.md](README.md)
* Increase version number in [configure.ac](configure.ac) in ```AC_INIT```
* [Increase Libtool version](http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html) number in [configure.ac](configure.ac) in ```AC_SUBST```
* Update to current year in [src/help.cpp](src/help.cpp)

### Release process

* Run tests using ```make check```
* Increase version number (see <a href="#versioning">Versioning</a>)
* Build statically linked primecout binaries and upload them to [https://bintray.com](https://bintray.com)
* Update _**Precompiled binaries**_ section in [README.md](README.md)
* Update [ChangeLog](ChangeLog)
* Tag the new release in git
* Create a new release tarball using ```make dist``` and upload it to [https://bintray.com](https://bintray.com)
