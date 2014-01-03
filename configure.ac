AC_INIT([primecount], [0.11], [kim.walisch@gmail.com])
AC_SUBST([primecount_lib_version], [0:0:0])

AC_PREREQ([2.62])
AC_LANG([C++])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Werror])
AC_PROG_CXX
AC_OPENMP
AC_CANONICAL_HOST
AM_PROG_AR
LT_INIT

# Force building shared library on Windows
# Reference: <http://article.gmane.org/gmane.comp.gnu.mingw.user/18727>
case "$host_os" in
    cygwin*) AM_CONDITIONAL(FORCE_SHARED_LIBRARY, true);;
    mingw*)  AM_CONDITIONAL(FORCE_SHARED_LIBRARY, true);;
    *)       AM_CONDITIONAL(FORCE_SHARED_LIBRARY, false);;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
