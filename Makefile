##############################################################################
# Makefile for primecount (console version & library)
#
# Author:          Kim Walisch
# Contact:         kim.walisch@gmail.com
# Created:         09 June 2013
# Last modified:   01 July 2013
#
# Project home:    https://github.com/kimwalisch/primecount
##############################################################################

TARGET   := primecount
CXX      := c++
CXXFLAGS := -O2
BINDIR   := bin
INCDIR   := include
LIBDIR   := lib
OBJDIR   := obj
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
  phi.o

PRIMECOUNT_HEADERS := \
  $(INCDIR)/primecount.h \
  $(wildcard \
    $(SRCDIR)/*.h)

#-----------------------------------------------------------------------------
# Needed to suppress output while checking system features
#-----------------------------------------------------------------------------

NO_STDOUT := 1> /dev/null
NO_STDERR := 2> /dev/null
NO_OUTPUT := $(NO_STDOUT) $(NO_STDERR)

#-----------------------------------------------------------------------------
# Find the compiler's OpenMP flag
#-----------------------------------------------------------------------------

OPENMP_PROGRAM := '\#include <omp.h>\n int main() { return _OPENMP; }'

is-openmp = $(shell command -v $(CXX) $(NO_OUTPUT) && \
                    printf $(OPENMP_PROGRAM) | \
                    $(CXX) $(CXXFLAGS) $1 -xc++ -c -o /dev/null - $(NO_STDERR) && \
                    echo successfully compiled!)

ifeq ($(call is-openmp,),)
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
# make            -> libprimecount.a
# make SHARED=yes -> libprimecount.(so|dylib)
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
# Default targets
#-----------------------------------------------------------------------------

.PHONY: all

all: bin lib

#-----------------------------------------------------------------------------
# Create and clean output directories
#-----------------------------------------------------------------------------

.PHONY: make_dir clean

make_dir:
	@mkdir -p $(BINDIR) $(LIBDIR) $(OBJDIR)

clean:
	rm -rf $(BINDIR) $(LIBDIR) $(OBJDIR)

#-----------------------------------------------------------------------------
# Compilation rules
#-----------------------------------------------------------------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) $(FPIC) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# Build the command-line program
#-----------------------------------------------------------------------------

.PHONY: bin bin_obj

bin: make_dir bin_obj

bin_obj: $(addprefix $(OBJDIR)/, $(PRIMECOUNT_OBJECTS))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$(TARGET) $^ -lprimesieve

#-----------------------------------------------------------------------------
# Build libprimecount
#-----------------------------------------------------------------------------

LIB_OBJECTS = $(addprefix \
                $(if $(FPIC),$(LIBDIR)/,$(OBJDIR)/), \
                  $(LIBPRIMECOUNT_OBJECTS))

.PHONY: lib lib_obj

lib: make_dir lib_obj

lib_obj: $(LIB_OBJECTS)
ifneq ($(SHARED),)
	$(CXX) $(strip $(CXXFLAGS) $(FPIC) $(SOFLAG)) -o $(LIBDIR)/$(LIBRARY) $^
else
	ar rcs $(LIBDIR)/$(LIBRARY) $^
endif

#-----------------------------------------------------------------------------
# `make check` runs correctness tests
#-----------------------------------------------------------------------------

.PHONY: check test

check test: bin
	$(BINDIR)/./$(TARGET) --test

#-----------------------------------------------------------------------------
# Install & uninstall targets
#-----------------------------------------------------------------------------

.PHONY: install uninstall

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
	@echo "make SHARED=yes                          Build a shared libprimecount library (using c++)"
	@echo "sudo make install                        Install primecount and libprimecount to /usr/local or /usr"
	@echo "sudo make install PREFIX=/path           Specify a custom installation path"
	@echo "sudo make uninstall                      Completely remove primecount and libprimecount"
	@echo "make help                                Print this help menu"
