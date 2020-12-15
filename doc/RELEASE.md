# How to release a new primecount version

* Update the [ChangeLog](../ChangeLog) file.
* Run ```scripts/update_version.sh``` in the root directory of the primecount git repo to update the version number.
* Run ```scripts/update_libprimesieve.sh``` in the root directory of the primecount git repo to update to the latest libprimesieve version.
* Run ```scripts/build_mingw64.sh``` (on Windows using MSYS2/MinGW-w64) in the root directory of the primecount git repo to build primecount release binaries for Windows (readily packaged as ```.zip``` archives).
* Go to the [GitHub website and do the release](https://github.com/kimwalisch/primecount/releases). The release title should be primecount-X.Y and the tag name should be vX.Y (e.g. primecount-1.0 and v1.0). For primecount-backup use e.g. primecount-backup-1.0 and v1.0-backup.
* First you should do the release of primecount-backup (```backup3``` branch) and afterwards the release of primecount (```master``` branch). In order to avoid trouble with the ordering of these 2 releases make sure that the ```master``` branch has the most recent commit e.g. by updating the date in the ```ChangeLog``` file on the master branch only.
