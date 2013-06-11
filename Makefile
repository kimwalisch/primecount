##############################################################################
# Makefile for primecount (console version & library)
#
# Author:          Kim Walisch
# Contact:         kim.walisch@gmail.com
# Created:         9 June 2013
# Last modified:   11 June 2013
#
# Project home:    https://github.com/kimwalisch/primecount
##############################################################################

CXX      := c++
CXXFLAGS := -O2
TARGET   := primecount
BINDIR   := bin
LIBDIR   := lib
INCDIR   := include
LEGENDRE := src/legendre
MEISSEL  := src/meissel
PK       := src/Pk
PROGRAM  := src/program
TEST     := src/test
UTILS    := src/utils

PRIMECOUNT_HEADERS := \
  $(INCDIR)/primecount.h \
  $(UTILS)/ExpressionParser.h \
  $(UTILS)/isqrt.h \
  $(UTILS)/Next_N_Primes_Vector.h \
  $(UTILS)/PrimeSieveVector.h

LIBPRIMECOUNT_OBJECTS := \
  pi_meissel.o \
  pi_legendre.o \
  phi.o \
  P2.o

PRIMECOUNT_OBJECTS := \
  primecount.o \
  pi_meissel.o \
  pi_legendre.o \
  phi.o \
  P2.o

TEST_OBJECTS := \
  primecount_test.o \
  pi_meissel.o \
  pi_legendre.o \
  phi.o \
  P2.o

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
# By default build primecount (command-line program) and libprimecount
#-----------------------------------------------------------------------------

.PHONY: all

all: bin lib

#-----------------------------------------------------------------------------
# Build the primecount console application
#-----------------------------------------------------------------------------

BIN_PRIMECOUNT_OBJECTS := \
  $(addprefix $(BINDIR)/, $(PRIMECOUNT_OBJECTS))

BIN_TEST_OBJECTS := \
  $(addprefix $(BINDIR)/, $(TEST_OBJECTS))

.PHONY: bin bin_dir primecount primecount_test

bin: bin_dir primecount primecount_test

bin_dir:
	@mkdir -p $(BINDIR)

primecount: $(BIN_PRIMECOUNT_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$(TARGET) $^ -lprimesieve

primecount_test: $(BIN_TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$(TARGET)_test $^ -lprimesieve

$(BINDIR)/%.o: $(LEGENDRE)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BINDIR)/%.o: $(MEISSEL)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BINDIR)/%.o: $(PK)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BINDIR)/%.o: $(PROGRAM)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BINDIR)/%.o: $(TEST)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# Build libprimecount
#-----------------------------------------------------------------------------

LIB_CXXFLAGS := $(strip $(CXXFLAGS) $(FPIC))
LIB_OBJECTS  := \
  $(addprefix $(LIBDIR)/, \
    $(notdir $(LIBPRIMECOUNT_OBJECTS)))

.PHONY: lib lib_dir lib_obj

lib: lib_dir lib_obj

lib_dir:
	@mkdir -p $(LIBDIR)

lib_obj: $(LIB_OBJECTS)
ifneq ($(SHARED),)
	$(CXX) $(LIB_CXXFLAGS) $(SOFLAG) -o $(LIBDIR)/$(LIBRARY) $^
else
	ar rcs $(LIBDIR)/$(LIBRARY) $^
endif

$(LIBDIR)/%.o: $(LEGENDRE)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(LIBDIR)/%.o: $(MEISSEL)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(LIBDIR)/%.o: $(PK)/%.cpp $(PRIMECOUNT_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

#-----------------------------------------------------------------------------
# `make check` runs correctness tests
#-----------------------------------------------------------------------------

.PHONY: check test

check test: bin_dir primecount_test
	$(BINDIR)/./primecount_test

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
