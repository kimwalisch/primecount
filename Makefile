##############################################################################
# Makefile for primecount (console version & library)
#
# Author:          Kim Walisch
# Contact:         kim.walisch@gmail.com
# Created:         9 June 2013
# Last modified:   9 June 2013
#
# Project home:    https://github.com/kimwalisch/primecount
##############################################################################

CXX      := c++
CXXFLAGS := -O2
TARGET   := primecount
BINDIR   := bin
LIBDIR   := lib
LEGENDRE := src/legendre
MEISSEL  := src/meissel
PROGRAMS := src/programs
TEST     := src/test
UTILS    := src/utils

HEADERS := \
  $(wildcard $(LEGENDRE)/*.h) \
  $(wildcard $(MEISSEL)/*.h) \
  $(wildcard $(PROGRAMS)/*.h) \
  $(wildcard $(TEST)/*.h) \
  $(wildcard $(UTILS)/*.h)

LEGENDRE_OBJS := $(subst .cpp,.o, $(wildcard $(LEGENDRE)/*.cpp))
MEISSEL_OBJS  := $(subst .cpp,.o, $(wildcard  $(MEISSEL)/*.cpp))

#-----------------------------------------------------------------------------
# Needed to suppress output while checking system features
#-----------------------------------------------------------------------------

NO_STDOUT := 1> /dev/null
NO_STDERR := 2> /dev/null
NO_OUTPUT := $(NO_STDOUT) $(NO_STDERR)

#-----------------------------------------------------------------------------
# Find the compiler's OpenMP flag
#-----------------------------------------------------------------------------

is-openmp = $(shell command -v $(CXX) $(NO_OUTPUT) && \
                    echo 'int main() { return _OPENMP; }' | \
                    $(CXX) $(CXXFLAGS) $1 -xc++ -c -o /dev/null - $(NO_STDERR) && \
                    echo successfully compiled!)

ifeq ($(call is-openmp),)
  ifneq ($(call is-openmp,-openmp),)
    CXXFLAGS += -openmp
  else
    ifneq ($(call is-openmp,-fopenmp),)
      CXXFLAGS += -fopenmp
    endif
  endif
endif

#-----------------------------------------------------------------------------
# `make`            -> libprimecount.a
# `make SHARED=yes` -> libprimecount.(so|dylib)
#-----------------------------------------------------------------------------

ifeq ($(SHARED),)
  LIBRARY := lib$(TARGET).a
else
  ifneq ($(shell uname | grep -i darwin),)
    SOFLAG := -dynamiclib
    LIBRARY := lib$(TARGET).dylib
  else
    SOFLAG := -shared
    LIBRARY := lib$(TARGET).so
    ifeq ($(shell uname | egrep -i 'mingw|cygwin'),)
      FPIC := -fPIC
    endif
  endif
endif

#-----------------------------------------------------------------------------
# Build the primesieve console application
#-----------------------------------------------------------------------------

.PHONY: bin bin_dir pi_legendre pi_meissel pi_test

bin: bin_dir pi_legendre pi_meissel pi_test

bin_dir:
	@mkdir -p $(BINDIR)

pi_legendre: $(PROGRAMS)/pi_legendre.o $(LEGENDRE_OBJS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/pi_legendre $^ -lprimesieve

pi_meissel: $(PROGRAMS)/pi_meissel.o $(MEISSEL_OBJS) $(LEGENDRE_OBJS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/pi_meissel $^ -lprimesieve

pi_test: $(TEST)/pi_test.o $(MEISSEL_OBJS) $(LEGENDRE_OBJS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/pi_test $^ -lprimesieve

$(LEGENDRE)/%.o: $(LEGENDRE)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MEISSEL)/%.o: $(MEISSEL)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PROGRAMS)/%.o: $(PROGRAMS)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST)/%.o: $(TEST)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#-----------------------------------------------------------------------------
# `make check` runs correctness tests
#-----------------------------------------------------------------------------

.PHONY: check test

check test: bin_dir pi_test
	$(BINDIR)/./pi_test

#-----------------------------------------------------------------------------
# Common targets (all, clean, install, uninstall)
#-----------------------------------------------------------------------------

.PHONY: all clean install uninstall

all: bin lib

clean:
	rm -rf $(BINDIR) $(LIBDIR) $(shell find . -name '*.o*')

# requires sudo privileges
install:
ifneq ($(wildcard $(BINDIR)/pi_*),)
	@mkdir -p $(PREFIX)/bin
	cp -f $(BINDIR)/pi_* $(PREFIX)/bin
endif
ifneq ($(wildcard $(LIBDIR)/lib$(TARGET).*),)
	@mkdir -p $(PREFIX)/lib
	cp -f $(wildcard $(LIBDIR)/lib$(TARGET).*) $(PREFIX)/lib
	cp -Rf $(INCDIR) $(PREFIX)
  ifneq ($(wildcard $(LIBDIR)/lib$(TARGET).so),)
    ifneq ($(shell command -v ldconfig $(NO_STDERR)),)
		ldconfig $(PREFIX)/lib
    endif
  endif
endif

# requires sudo privileges
uninstall:
ifneq ($(wildcard $(PREFIX)/bin/pi_*),)
	rm -f $(PREFIX)/bin/pi_legendre
	@rm -f $(PREFIX)/bin/pi_legendre.exe
	rm -f $(PREFIX)/bin/pi_meissel
	@rm -f $(PREFIX)/bin/pi_meissel.exe
endif
ifneq ($(wildcard $(PREFIX)/include/$(TARGET)),)
	rm -rf $(PREFIX)/include/$(TARGET)
endif
ifneq ($(wildcard $(PREFIX)/lib/lib$(TARGET).*),)
  ifneq ($(wildcard $(PREFIX)/lib/lib$(TARGET).so),)
		rm -f $(wildcard $(PREFIX)/lib/lib$(TARGET).so)
    ifneq ($(shell command -v ldconfig $(NO_STDERR)),)
		ldconfig $(PREFIX)/lib
    endif
  else
	rm -f $(wildcard $(PREFIX)/lib/lib$(TARGET).*)
  endif
endif

