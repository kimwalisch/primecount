AC_INIT([primecount], [1.0], [kim.walisch@gmail.com])
AC_SUBST([primecount_lib_version], [1:0:0])

AC_PREREQ([2.62])
AC_LANG([C++])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_PROG_CXX
AC_CANONICAL_HOST
AC_OPENMP
AM_PROG_AR
LT_INIT

# Enable POPCNT instruction by default
AC_ARG_ENABLE([popcnt],
    AC_HELP_STRING([--enable-popcnt], [enable POPCNT instruction (default yes)]))

AS_IF([test "x$enable_popcnt" != "xno"], [
    AX_GCC_BUILTIN([__builtin_popcount])
    AX_GCC_BUILTIN([__builtin_popcountll])
    AS_IF([test "x$ax_cv_have___builtin_popcount" = "xyes"], [
        AX_POPCNT
    ])
])

# Disable assert macro by default
AC_ARG_ENABLE([assert],
    AC_HELP_STRING([--enable-assert], [enable assert macro (default no)]),
    [], [AC_DEFINE(NDEBUG, 1, [disable assertions])])

# Force building shared library on Windows
# Reference: <http://article.gmane.org/gmane.comp.gnu.mingw.user/18727>
case "$host_os" in
    cygwin*) AM_CONDITIONAL(FORCE_SHARED_LIBRARY, true);;
    mingw*)  AM_CONDITIONAL(FORCE_SHARED_LIBRARY, true);;
    *)       AM_CONDITIONAL(FORCE_SHARED_LIBRARY, false);;
esac

AC_SEARCH_LIBS([primesieve_generate_primes], [primesieve], [],
    [AC_MSG_ERROR([libprimesieve is missing, get it from <http://primesieve.org>])]
)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
