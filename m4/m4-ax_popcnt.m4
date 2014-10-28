# ===========================================================================
#   https://github.com/kimwalisch/primecount/blob/master/m4/m4-ax_popcnt.m4
# ===========================================================================
#
# SYNOPSIS
#
#   AX_POPCNT
#
# DESCRIPTION
#
#   Find if the CPU supports the POPCNT instruction by requesting cpuid. If
#   so then add -mpopcnt (or -msse4.2) to POPCNT_FLAGS if the compiler
#   supports it.
#
# ACKNOWLEDGEMENTS
#
#   The m4-ax_popcnt.m4 macro is a modified version of the m4-ax_ext.m4
#   macro written by Christophe Tournayre and Michael Petch.
#
#   This macro calls:
#
#     AC_SUBST(POPCNT_FLAGS)
#
#   And defines:
#
#     HAVE_POPCNT
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013 Michael Petch <mpetch@capp-sysware.com>
#   Copyright (c) 2014 Kim Walisch <kim.walisch@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 13

AC_DEFUN([AX_POPCNT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  case $host_cpu in
    i[[3456]]86*|x86_64*|amd64*)

      AC_REQUIRE([AX_GCC_X86_CPUID])

      AX_GCC_X86_CPUID(0x00000001)
      ecx=0
      if test "$ax_cv_gcc_x86_cpuid_0x00000001" != "unknown";
      then
        ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
      fi

      AC_CACHE_CHECK([whether POPCNT is supported], [ax_cv_have_popcnt_ext],
      [
        ax_cv_have_popcnt_ext=no
        if test "$((0x$ecx>>23&0x01))" = 1; then
          ax_cv_have_popcnt_ext=yes
        fi
      ])

      if test "$ax_cv_have_popcnt_ext" = no; then
        AC_MSG_NOTICE([NOTICE: Your CPU does not support the POPCNT instruction.])
      else
        AX_CHECK_COMPILE_FLAG(-mpopcnt, ax_cv_support_popcnt_ext=yes, [])
        if test x"$ax_cv_support_popcnt_ext" = x"yes"; then
          POPCNT_FLAGS="-mpopcnt -DHAVE_POPCNT=1"
        else
          AC_CACHE_CHECK([whether sse4.2 is supported], [ax_cv_have_sse42_ext],
          [
            ax_cv_have_sse42_ext=no
            if test "$((0x$ecx>>20&0x01))" = 1; then
              ax_cv_have_sse42_ext=yes
            fi
          ])
          if test "$ax_cv_have_sse42_ext" = yes; then
            AX_CHECK_COMPILE_FLAG(-msse4.2, ax_cv_support_sse42_ext=yes, [])
            if test x"$ax_cv_support_sse42_ext" = x"yes"; then
              POPCNT_FLAGS="-msse4.2 -DHAVE_POPCNT=1"
            fi
          fi
          if test x"$POPCNT_FLAGS" = x; then
            AC_MSG_NOTICE([NOTICE: Your compiler does not support the POPCNT instruction.])
          fi
        fi
      fi
    ;;

    *) # Non x86 CPU architectures
      for i in $($CXX -march=native -E -v - </dev/null 2>&1 | grep -o " -\w*popc\w* " || $CXX --help=target 2>&1 | grep -o " -\w*popc\w* " || $CXX --target-help 2>&1 | grep -o " -\w*popc\w* "); do
        POPCNT_FLAG="$i"
        if test x"$POPCNT_FLAG" != x; then
          AX_CHECK_COMPILE_FLAG($POPCNT_FLAG, ax_cv_support_popcnt_ext=yes, [])
          if test x"$ax_cv_support_popcnt_ext" = x"yes"; then
            POPCNT_FLAGS="$POPCNT_FLAGS $POPCNT_FLAG"
          fi
        fi
      done
      if test x"$POPCNT_FLAGS" != x; then
        POPCNT_FLAGS="$POPCNT_FLAGS -DHAVE_POPCNT=1"
      else
        AC_MSG_NOTICE([NOTICE: Your compiler and/or CPU do not support the POPCNT instruction.])
      fi
    ;;
  esac

  AC_SUBST(POPCNT_FLAGS)
])
