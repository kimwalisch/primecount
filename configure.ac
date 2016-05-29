AC_INIT([primecount], [3.3], [kim.walisch@gmail.com])
AC_SUBST([primecount_lib_version], [4:6:0])

AC_PREREQ([2.62])
AC_LANG([C++])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_PROG_CXX
AC_CANONICAL_HOST
AC_SUBST(AR_FLAGS, [cr])
AC_OPENMP
AM_PROG_AR
LT_INIT

# Check if OS is Windows
AM_CONDITIONAL(IS_WINDOWS, test "$(uname 2>/dev/null | egrep -i 'windows|cygwin|mingw|msys')" != "")

# Disable assert macro by default
AC_ARG_ENABLE([assert],
    AC_HELP_STRING([--enable-assert],
                   [enable assert macro (default no)]),
    [], [AC_DEFINE(NDEBUG, 1, [disable assertions])])

# Disable MPI by default
AC_ARG_ENABLE([mpi],
    AC_HELP_STRING([--enable-mpi],
                   [enable MPI (Message Passing Interface)]),
    [], [])

AS_IF([test "x$enable_mpi" = "xyes"], [
    AC_CHECK_PROGS(MPICXX, mpic++,)
    test -z "$MPICXX" && AC_MSG_ERROR([couldn't find MPI library])
    CXX="$MPICXX"
    AC_DEFINE([HAVE_MPI], 1, [Use MPI])
])

AM_CONDITIONAL(HAVE_MPI, test x$enable_mpi = xyes)

# Enable POPCNT instruction by default
AC_ARG_ENABLE([popcnt],
    AC_HELP_STRING([--enable-popcnt],
                   [enable POPCNT instruction (default yes)]))

AS_IF([test "x$enable_popcnt" != "xno"], [
    AC_MSG_CHECKING(for __builtin_popcount)
    AC_TRY_LINK(, [int x = __builtin_popcount(123);],
                  [AC_MSG_RESULT(yes); builtin_popcount=yes; AC_DEFINE(HAVE___BUILTIN_POPCOUNT)],
                  [AC_MSG_RESULT(no)])

    AC_MSG_CHECKING(for __builtin_popcountll)
    AC_TRY_LINK(, [long long x = __builtin_popcountll(123);],
                  [AC_MSG_RESULT(yes); builtin_popcount=yes; AC_DEFINE(HAVE___BUILTIN_POPCOUNTLL)],
                  [AC_MSG_RESULT(no)])

    AS_IF([test "x$builtin_popcount" = "xyes"], [
        AX_POPCNT
    ])
])

AC_MSG_CHECKING(for int128_t)
AC_TRY_LINK([#include <stdint.h>], [int128_t x = 0;],
            [AC_MSG_RESULT(yes); int128=yes; AC_DEFINE(HAVE_INT128_T)],
            [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __int128_t)
AC_TRY_LINK(, [__int128_t x = 0;],
              [AC_MSG_RESULT(yes); int128=yes; AC_DEFINE(HAVE___INT128_T)],
              [AC_MSG_RESULT(no)])

AM_CONDITIONAL(HAVE_INT128_T, test x$int128 = xyes)

# Enable libdivide by default
AC_ARG_ENABLE([libdivide],
    AC_HELP_STRING([--enable-libdivide],
                   [enable libdivide (optimized integer division)]))

AS_IF([test "x$enable_libdivide" != "xno"], [
    AC_MSG_CHECKING(for libdivide)
    AC_TRY_RUN([
    #include "include/libdivide.h"
    #include <stdint.h>
    #include <limits>

    int main()
    {
        uint64_t tests[20] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            63, 101, 511,
            1 << 5, 1 << 9, 1 << 20,
            std::numeric_limits<int32_t>::max(),
            std::numeric_limits<uint32_t>::max(),
            std::numeric_limits<int64_t>::max(),
            std::numeric_limits<uint64_t>::max() };

        for (int i = 0; i < 20; i++)
        {
            uint64_t x = tests[i];

            for (int j = 2; j < 1000; j++)
            {
                libdivide::divider<uint64_t, libdivide::BRANCHFREE> divider(j);
                if (x / j != x / divider)
                    return 1;
            }
        }

        return 0;
    }
    ], [AC_MSG_RESULT(yes); libdivide=yes; AC_DEFINE(HAVE_LIBDIVIDE)],
       [AC_MSG_RESULT(no)])
])

AM_CONDITIONAL(HAVE_LIBDIVIDE, test x$libdivide = xyes)

AC_SEARCH_LIBS([primesieve_generate_primes], [primesieve], [],
    [AC_MSG_ERROR([libprimesieve is missing, get it from <http://primesieve.org>])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
