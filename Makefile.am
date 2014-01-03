ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libprimecount.la
libprimecount_la_CXXFLAGS = $(OPENMP_CXXFLAGS)
libprimecount_la_LDFLAGS = -version-info @primecount_lib_version@

if FORCE_SHARED_LIBRARY
libprimecount_la_LDFLAGS += -no-undefined 
endif

bin_PROGRAMS = primecount
primecount_LDADD = libprimecount.la -lprimesieve

# If installation directory is /usr/local/lib
# execute ldconfig after installation
install-exec-hook:
	test "$(prefix)" = "/usr/local" && command -v ldconfig 2>/dev/null && ldconfig /usr/local/lib 2>/dev/null || true

check: all
	./primecount --test

include_HEADERS = \
	include/primecount.hpp

EXTRA_DIST = \
	images \
	README.md

libprimecount_la_SOURCES = \
  src/Li.cpp \
  src/nth_prime.cpp \
  src/Pk.cpp \
  src/pi_primesieve.cpp \
  src/pi_meissel.cpp \
  src/pi_legendre.cpp \
  src/pi_lehmer.cpp \
  src/pi_lmo_simple.cpp \
  src/pi_lmo.cpp \
  src/pi.cpp \
  src/phi.cpp \
  src/PhiTiny.cpp

primecount_SOURCES = \
  src/cmdoptions.cpp \
  src/help.cpp \
  src/primecount.cpp \
  src/test.cpp
