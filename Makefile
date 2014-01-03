##############################################################################
# Makefile for primecount (console version & library)
# Project home:    https://github.com/kimwalisch/primecount
##############################################################################

TARGET   := primecount
CXX      := c++
CXXFLAGS := -O2 -DNDEBUG
BINDIR   := bin
INCDIR   := include
LIBDIR   := lib
OBJDIR   := obj
FPICDIR  := obj/fpic
SRCDIR   := src

LIB_OBJECTS := \
  $(OBJDIR)/Li.o \
  $(OBJDIR)/nth_prime.o \
  $(OBJDIR)/Pk.o \
  $(OBJDIR)/pi_primesieve.o \
  $(OBJDIR)/pi_meissel.o \
  $(OBJDIR)/pi_legendre.o \
  $(OBJDIR)/pi_lehmer.o \
  $(OBJDIR)/pi_lmo_simple.o \
  $(OBJDIR)/pi_lmo.o \
  $(OBJDIR)/pi.o \
  $(OBJDIR)/phi.o \
  $(OBJDIR)/PhiTiny.o

BIN_OBJECTS := \
  $(LIB_OBJECTS) \
  $(OBJDIR)/cmdoptions.o \
  $(OBJDIR)/help.o \
  $(OBJDIR)/primecount.o \
  $(OBJDIR)/test.o

HEADERS := \
  $(INCDIR)/primecount.hpp \
  $(wildcard $(SRCDIR)/*.hpp)

#-----------------------------------------------------------------------------
# Needed to suppress output while checking system features
#-----------------------------------------------------------------------------

NO_STDOUT := 1> /dev/null
NO_STDERR := 2> /dev/null
NO_OUTPUT := $(NO_STDOUT) $(NO_STDERR)

#-----------------------------------------------------------------------------
# Find the compiler's OpenMP flag
#-----------------------------------------------------------------------------

ifneq ($(OPENMP),no)
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
# make        -> libprimesieve.a
# make shared -> libprimesieve.(so|dylib)
#-----------------------------------------------------------------------------

ifneq ($(shell uname | grep -i darwin),)
  SOFLAG := -dynamiclib
  SHARED_LIBRARY := lib$(TARGET).dylib
else
  SOFLAG := -shared
  SHARED_LIBRARY := lib$(TARGET).so
  ifeq ($(shell uname | egrep -i 'mingw|cygwin'),)
    FPIC := -fPIC
  endif
endif

#-----------------------------------------------------------------------------
# Default targets
#-----------------------------------------------------------------------------

.PHONY: default all

default: bin static

all: bin lib

#-----------------------------------------------------------------------------
# Create and clean output directories
#-----------------------------------------------------------------------------

.PHONY: make_dir clean

make_dir:
	@mkdir -p $(BINDIR) $(LIBDIR) $(OBJDIR) $(FPICDIR)

clean:
	rm -rf $(BINDIR) $(LIBDIR) $(OBJDIR) $(FPICDIR)

#-----------------------------------------------------------------------------
# Compilation rules
#-----------------------------------------------------------------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(FPICDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(FPIC) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# Build the command-line program
#-----------------------------------------------------------------------------

.PHONY: bin bin_obj

bin: make_dir bin_obj

bin_obj: $(BIN_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$(TARGET) $^ -lprimesieve

#-----------------------------------------------------------------------------
# Build libprimecount
#-----------------------------------------------------------------------------

SHARED_OBJECTS = $(if $(FPIC), \
                   $(subst $(OBJDIR),$(FPICDIR),$(LIB_OBJECTS)), \
                     $(LIB_OBJECTS))

.PHONY: lib static shared static_obj shared_obj

lib: static shared

static: make_dir static_obj
shared: make_dir shared_obj

static_obj: $(LIB_OBJECTS)
	$(AR) rcs $(LIBDIR)/lib$(TARGET).a $^

shared_obj: $(SHARED_OBJECTS)
	$(CXX) $(strip $(CXXFLAGS) $(FPIC) $(SOFLAG)) -o $(LIBDIR)/$(SHARED_LIBRARY) $^ -lprimesieve

#-----------------------------------------------------------------------------
# Run integration tests
#-----------------------------------------------------------------------------

.PHONY: check test

check test: bin
	$(BINDIR)/./$(TARGET) --test

#-----------------------------------------------------------------------------
# Install & uninstall targets
#-----------------------------------------------------------------------------

.PHONY: install uninstall

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
      ifeq ($(firstword $(subst /, ,$(PREFIX))),usr)
			ldconfig $(PREFIX)/lib
      endif
    endif
  endif
endif

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
		rm -f $(wildcard $(PREFIX)/lib/lib$(TARGET).*)
    ifneq ($(shell command -v ldconfig $(NO_STDERR)),)
      ifeq ($(firstword $(subst /, ,$(PREFIX))),usr)
			ldconfig $(PREFIX)/lib
      endif
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
	@echo "make                                    Build primecount and static libprimecount"
	@echo "make shared                             Build shared libprimecount, requires shared libprimesieve"
	@echo "make CXX=icpc CXXFLAGS=\"-O2 -openmp\"    Specify a custom C++ compiler, here icpc"
	@echo "make check                              Run integration tests"
	@echo "sudo make install                       Install primecount and libprimecount to /usr[/local]"
	@echo "sudo make install PREFIX=/path          Specify a custom installation path"
	@echo "sudo make uninstall                     Completely remove primecount and libprimecount"
	@echo "make clean                              Clean the output directories (bin, lib, ...)"
	@echo "make help                               Print this help menu"
