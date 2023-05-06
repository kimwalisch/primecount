# How to release a new primecount version

* Update the [ChangeLog](../ChangeLog) file.
* Run ```scripts/update_version.sh``` in the root directory of the primecount git repo to update the version number. This script takes the new version number as a parameter e.g.: ```scripts/./update_version.sh 1.2```
* Run ```scripts/update_libprimesieve.sh``` in the root directory of the primecount git repo to update to the latest libprimesieve version.
* Go to the [GitHub website and do the release](https://github.com/kimwalisch/primecount/releases). The release title should be primecount-X.Y and the tag name should be vX.Y (e.g. primecount-1.0 and v1.0).

# Optionally

* Run ```scripts/build_mingw64_x64.sh``` (on Windows using MSYS2/MinGW-w64) in the root directory of the primecount git repo to build primecount release binaries for Windows (readily packaged as ```.zip``` archives).
* When doing the release on the GitHub website, add the new primecount-win64 zip archive as an artifact to the release.
