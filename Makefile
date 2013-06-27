##############################################################################
# Makefile for primecount (console version & library)
#
# Author:          Kim Walisch
# Contact:         kim.walisch@gmail.com
# Created:         9 June 2013
# Last modified:   27 June 2013
#
# Project home:    https://github.com/kimwalisch/primecount
##############################################################################

CXX      := c++
CXXFLAGS := -O2
TARGET   := primecount
BINDIR   := bin
LIBDIR   := lib
INCDIR   := include
SRCDIR   := src

PRIMECOUNT_OBJECTS := \
  cmdoptions.o \
  help.o \
  Li.o \
  nth_prime.o \
  primecount.o \
  pi_primesieve.o \
  pi_meissel.o \
  pi_legendre.o \
  pi_lehmer.o \
  pi.o \
  phi.o \
  test.o

LIBPRIMECOUNT_OBJECTS := \
  Li.o \
  nth_prime.o \
  pi_primesieve.o \
  pi_meissel.o \
  pi_legendre.o \
  pi_lehmer.o \
  pi.o \
  phi.o \
  test.o

PRIMECOUNT_HEADERS := \
  $(INCDIR)/primecount.h \
  $(SRCDIR)/phi.h \
  $(SRCDIR)/pi_bsearch.h \
  $(SRCDIR)/cmdoptions.h \
  $(SRCDIR)/ExpressionParser.h \
  $(SRCDIR)/isqrt.h

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
# Default installation path
#-----------------------------------------------------------------------------

PREFIX := /usr

ifneq ($(shell uname | grep -i linux),)
  PREFIX := /usr/local
endif
ifneq ($(shell uname | grep -i mingw),)
  PREFIX := /mingw
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
# By default build command-line programs and libprimecount
#-----------------------------------------------------------------------------

.PHONY: all

all: bin lib

#-----------------------------------------------------------------------------
# Build the command-line programs
#-----------------------------------------------------------------------------

.PHONY: bin bin_dir bin_obj

bin: bin_dir bin_obj

bin_dir:
	@mkdir -p $(BINDIR)

bin_obj: $(addprefix $(BINDIR)/, $(PRIMECOUNT_OBJECTS))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$(TARGET) $^ -lprimesieve

$(BINDIR)/%.o: $(SRCDIR)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# Build libprimecount
#-----------------------------------------------------------------------------

LIB_CXXFLAGS := $(strip $(CXXFLAGS) $(FPIC))

.PHONY: lib lib_dir lib_obj

lib: lib_dir lib_obj

lib_dir:
	@mkdir -p $(LIBDIR)

lib_obj: $(addprefix $(LIBDIR)/, $(LIBPRIMECOUNT_OBJECTS))
ifneq ($(SHARED),)
	$(CXX) $(LIB_CXXFLAGS) $(SOFLAG) -o $(LIBDIR)/$(LIBRARY) $^
else
	ar rcs $(LIBDIR)/$(LIBRARY) $^
endif

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# `make check` runs correctness tests
#-----------------------------------------------------------------------------

.PHONY: check test

check test: bin
	$(BINDIR)/./$(TARGET) --test

#-----------------------------------------------------------------------------
# Common targets (clean, install, uninstall)
#-----------------------------------------------------------------------------

.PHONY: clean install uninstall

clean:
	rm -rf $(BINDIR) $(LIBDIR)

# requires sudo privileges
install:
ifneq ($(wildcard $(BINDIR)/$(TARGET)*),)
	@mkdir -p $(PREFIX)/bin
	cp -f $(BINDIR)/$(TARGET) $(PREFIX)/bin
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
ifneq ($(wildcard $(PREFIX)/bin/$(TARGET)*),)
	rm -f $(PREFIX)/bin/$(TARGET)
	@rm -f $(PREFIX)/bin/$(TARGET).exe
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

#-----------------------------------------------------------------------------
# Makefile help menu
#-----------------------------------------------------------------------------

.PHONY: help

help:
	@echo ----------------------------------------------
	@echo ---------- primecount build options ----------
	@echo ----------------------------------------------
	@echo "make                                     Build the primecount console application (using c++)"
	@echo "make CXX=icpc CXXFLAGS=\"-O2 -openmp\"     Specify a custom C++ compiler, here icpc"
	@echo "make check                               Test primecount for correctness"
	@echo "make clean                               Clean the output directories (bin, lib, ...)"
	@echo "make lib                                 Build a static libprimecount library (using c++)"
	@echo "make lib SHARED=yes                      Build a shared libprimecount library (using c++)"
	@echo "sudo make install                        Install primecount and libprimecount to /usr/local or /usr"
	@echo "sudo make install PREFIX=/path           Specify a custom installation path"
	@echo "sudo make uninstall                      Completely remove primecount and libprimecount"
	@echo "make help                                Print this help menu"
