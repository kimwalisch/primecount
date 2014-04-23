Hacking on primecount
=====================

Add a new prime counting implementation
---------------------------------------
* Add new function signature e.g. ```pi_lmo(x)``` to [include/primecount.hpp](include/primecount.hpp)
* Create source file ```src/pi_lmo.cpp``` with the new function
* Add ```src/pi_lmo.cpp``` to [Makefile.am](Makefile.am)
* Add ```src\pi_lmo.cpp``` to [Makefile.msvc](Makefile.msvc)
* Add ```OPTION_LMO``` to [src/cmdoptions.hpp](src/cmdoptions.hpp)
* Add ```OPTION_LMO``` to optionMap in [src/cmdoptions.cpp](src/cmdoptions.cpp)
* Add ```pi_lmo(x)``` to [src/primecount.cpp](src/primecount.cpp)
* Add ```pi_lmo(x)``` to ```bool test()``` in [src/test.cpp](src/test.cpp)
* Add ```--lmo``` command-line option summary to [src/help.cpp](src/help.cpp)

Versioning
----------

* Increase version number in [include/primecount.hpp](include/primecount.hpp)
* Increase version number in _**Build instructions**_ section in [README.md](README.md)
* Increase version number in [configure.ac](configure.ac) in ```AC_INIT```
* [Increase Libtool version](http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html) number in [configure.ac](configure.ac) in ```AC_SUBST```
* Update to current year in [src/help.cpp](src/help.cpp)

Release process
---------------

* Increase version number (see [Versioning](#versioning"))
* Build statically linked primecout binaries for Linux and Windows
* Upload binaries to [bintray.com](bintray.com)
* Update _**Precompiled binaries**_ section in [README.md](README.md)
* Update [ChangeLog](ChangeLog)
