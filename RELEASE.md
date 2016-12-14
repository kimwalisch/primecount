Release process
===============

Clone or fork primecount:
```sh
git clone git://github.com/kimwalisch/primecount.git
```

Build primecount using:
```sh
cmake .
make -j8
```

### Adding a new prime counting function implementation

* Add new function signature e.g. ```pi_lmo(x)``` to [include/primecount.hpp](include/primecount-internal.hpp)
* Create source file ```src/lmo/pi_lmo.cpp``` with the new function
* Add ```src/lmo/pi_lmo.cpp``` to [CMakeLists.txt](CMakeLists.txt)
* Add ```pi_lmo(x)``` to [src/app/main.cpp](src/app/main.cpp)
* Add ```pi_lmo(x)``` to ```bool test()``` in [src/test.cpp](src/test.cpp)
* Add ```--lmo``` command-line option summary to [src/app/help.cpp](src/app/help.cpp)
* Add ```OPTION_LMO``` to [src/app/cmdoptions.hpp](src/app/cmdoptions.hpp)
* Add ```OPTION_LMO``` to optionMap in [src/app/cmdoptions.cpp](src/app/cmdoptions.cpp)

### Versioning

Increase the primecount version before each new release. The
```update_version.sh``` script automatically updates the version
(and release date) in all source files.

```sh
cd scripts

# Usage example: update primecount version to 3.5
./update_version.sh 3.5
```
The libprimecount version must be updated manually in ```CMakeLists.txt```.

### Release process

* Run tests: ```./primecount --test```
* Increase version number (see <a href="#versioning">Versioning</a>)
* Build statically linked primecout binaries and upload them to [https://bintray.com/kimwalisch/primecount](https://bintray.com/kimwalisch/primecount)
* Update _**Precompiled binaries**_ section in [README.md](README.md)
* Update [ChangeLog](ChangeLog)
* Tag the new release in git
* Create new release tarball using ```make dist``` and upload it to [https://bintray.com/kimwalisch/primecount](https://bintray.com/kimwalisch/primecount)
