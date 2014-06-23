ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libprimecount.la
libprimecount_la_CXXFLAGS = $(OPENMP_CXXFLAGS)
libprimecount_la_LDFLAGS = -version-info @primecount_lib_version@

if FORCE_SHARED_LIBRARY
libprimecount_la_LDFLAGS += -no-undefined 
endif

bin_PROGRAMS = primecount
primecount_LDADD = libprimecount.la

# If installation directory is /usr/local/lib
# execute ldconfig after installation
install-exec-hook:
	test "$(prefix)" = "/usr/local" && command -v ldconfig 2>/dev/null && ldconfig /usr/local/lib 2>/dev/null || true

check: all
	./primecount --test

include_HEADERS = \
	include/primecount.hpp

EXTRA_DIST = \
	docs/computing-special-leaves.md \
	src/has_openmp.cpp \
	Makefile.msvc \
	README.md

libprimecount_la_SOURCES = \
	src/balance_S2_load.cpp \
	src/Li.cpp \
	src/nth_prime.cpp \
	src/PhiCache.cpp \
	src/P2.cpp \
	src/P3.cpp \
	src/pi_primesieve.cpp \
	src/pi_meissel.cpp \
	src/pi_legendre.cpp \
	src/pi_lehmer.cpp \
	src/phi.cpp \
	src/PhiTiny.cpp \
	src/pmath.cpp \
	src/popcount.cpp \
	src/primecount.cpp \
	src/test.cpp \
	src/deleglise-rivat/pi_deleglise_rivat1.cpp \
	src/deleglise-rivat/pi_deleglise_rivat2.cpp \
	src/deleglise-rivat/pi_deleglise_rivat3.cpp \
	src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp \
	src/deleglise-rivat/pi_deleglise_rivat_parallel2.cpp \
	src/lmo/FactorTable.cpp \
	src/lmo/init_square_free.cpp \
	src/lmo/pi_lmo1.cpp \
	src/lmo/pi_lmo2.cpp \
	src/lmo/pi_lmo3.cpp \
	src/lmo/pi_lmo4.cpp \
	src/lmo/pi_lmo5.cpp \
	src/lmo/pi_lmo6.cpp \
	src/lmo/pi_lmo_parallel1.cpp \
	src/lmo/pi_lmo_parallel2.cpp \
	src/lmo/pi_lmo_parallel3.cpp \
	src/lmo/pi_lmo_parallel4.cpp \
	src/lmo/S1.cpp \
	include/aligned_vector.hpp \
	include/balance_S2_load.hpp \
	include/bit_sieve.hpp \
	include/FactorTable.hpp \
	include/init_square_free.hpp \
	include/primecount-internal.hpp \
	include/pmath.hpp \
	include/PhiCache.hpp \
	include/PhiTiny.hpp \
	include/tos_counters.hpp \
	include/utils.hpp

primecount_SOURCES = \
	src/app/main.cpp \
	src/app/help.cpp \
	src/app/cmdoptions.cpp \
	src/app/cmdoptions.hpp \
	src/app/calculator.hpp
