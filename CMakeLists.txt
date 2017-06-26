cmake_minimum_required(VERSION 3.1)
project(primecount)
set(PRIMECOUNT_VERSION_MAJOR 3)
set(PRIMECOUNT_VERSION_MINOR 8)
set(PRIMECOUNT_VERSION "${PRIMECOUNT_VERSION_MAJOR}.${PRIMECOUNT_VERSION_MINOR}")
set(CMAKE_BUILD_TYPE Release)

# Build options ######################################################

option(WITH_POPCNT       "Enable x86 POPCNT instruction"  ON)
option(WITH_LIBDIVIDE    "Use libdivide.h"                ON)
option(WITH_OPENMP       "Enable OpenMP support"          ON)
option(WITH_MPI          "Enable MPI support"             OFF)
option(BUILD_PRIMECOUNT  "Build primecount binary"        ON)
option(BUILD_SHARED_LIBS "Build shared libprimecount"     OFF)
option(BUILD_STATIC_LIBS "Build static libprimecount"     ON)
option(BUILD_TESTS       "Build test programs"            OFF)

if(WIN32)
    set(BUILD_SHARED_LIBS OFF)
endif()

# Compiler must support C++11 or later ###############################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Check* #####################################################

include(CheckCXXCompilerFlag)
include(CheckTypeSize)

# primecount binary source files #####################################

set(BIN_SRC src/app/cmdoptions.cpp
            src/app/help.cpp
            src/app/main.cpp)

# primecount library source files ####################################

set(LIB_SRC src/backup.cpp
            src/BitSieve.cpp
            src/FactorTable.cpp
            src/Li.cpp
            src/P2.cpp
            src/P3.cpp
            src/PhiTiny.cpp
            src/PiTable.cpp
            src/S1.cpp
            src/LoadBalancer.cpp
            src/S2Status.cpp
            src/Wheel.cpp
            src/generate.cpp
            src/nth_prime.cpp
            src/phi.cpp
            src/pi_legendre.cpp
            src/pi_lehmer.cpp
            src/pi_meissel.cpp
            src/pi_primesieve.cpp
            src/primecount.cpp
            src/print.cpp
            src/test.cpp
            src/lmo/pi_lmo1.cpp
            src/lmo/pi_lmo2.cpp
            src/lmo/pi_lmo3.cpp
            src/lmo/pi_lmo4.cpp
            src/lmo/pi_lmo5.cpp
            src/lmo/pi_lmo_parallel.cpp
            src/deleglise-rivat/S2_hard.cpp
            src/deleglise-rivat/S2_trivial.cpp
            src/deleglise-rivat/pi_deleglise_rivat1.cpp
            src/deleglise-rivat/pi_deleglise_rivat2.cpp
            src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp)

# Check for MPI (Message Passing Interface) ##########################

if(WITH_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    add_definitions(-DHAVE_MPI)

    set(LIB_SRC ${LIB_SRC}
                src/mpi/P2_mpi.cpp
                src/mpi/S2_hard_mpi.cpp
                src/mpi/MpiLoadBalancer.cpp
                src/mpi/MpiMsg.cpp)
endif()

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Check if int128_t type exists ######################################

check_type_size("int128_t" INT128_T)
if(HAVE_INT128_T)
    add_definitions(-DHAVE_INT128_T)
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/pi_deleglise_rivat_parallel2.cpp)
else()
    check_type_size("__int128_t" __INT128_T)
    if(HAVE___INT128_T)
        add_definitions(-DHAVE___INT128_T)
        set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/pi_deleglise_rivat_parallel2.cpp)
    endif()
endif()

# Check if libdivide.h works #########################################

if(WITH_LIBDIVIDE)
    message(STATUS "Performing Test libdivide")

    try_run(RUN_LIBDIVIDE COMPILE_LIBDIVIDE
            ${CMAKE_BINARY_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/test/libdivide.cpp
            COMPILE_DEFINITIONS -I${CMAKE_CURRENT_SOURCE_DIR}/include)

    if(RUN_LIBDIVIDE EQUAL 0)
        message(STATUS "Performing Test libdivide - Success")
    else()
        message(STATUS "Performing Test libdivide - Failed")
    endif()
endif()

# Include S2_easy* source files ######################################

if(RUN_LIBDIVIDE EQUAL 0)
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
    if(WITH_MPI)
        set(LIB_SRC ${LIB_SRC} src/mpi/S2_easy_mpi_libdivide.cpp)
    endif()
else()
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy.cpp)
    if(WITH_MPI)
        set(LIB_SRC ${LIB_SRC} src/mpi/S2_easy_mpi.cpp)
    endif()
endif()

# Check -mpopcnt compiler flag #######################################

if(WITH_POPCNT)
    set(COPY_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS -Werror)
    check_cxx_compiler_flag(-mpopcnt mpopcnt)
    set(CMAKE_CXX_FLAGS "${COPY_CXX_FLAGS}")
    if(mpopcnt)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif()
endif()

# Testing ############################################################

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# libprimesieve ######################################################

option(BUILD_PRIMESIEVE "Build primesieve binary" OFF)
set(BUILD_TESTS OFF CACHE BOOL "Build primesieve tests" FORCE)

add_subdirectory(lib/primesieve)

# libprimecount ######################################################

add_library(libprimecount ${LIB_SRC})
set_target_properties(libprimecount PROPERTIES OUTPUT_NAME primecount)
target_link_libraries(libprimecount libprimesieve ${MPI_CXX_LIBRARIES})

target_include_directories(libprimecount PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if(BUILD_SHARED_LIBS)
    set_target_properties(libprimecount PROPERTIES SOVERSION ${PRIMECOUNT_VERSION_MAJOR})
    set_target_properties(libprimecount PROPERTIES VERSION ${PRIMECOUNT_VERSION})
endif()

install(TARGETS libprimecount
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# static_libprimecount ###############################################

if(BUILD_STATIC_LIBS AND BUILD_SHARED_LIBS)
    add_library(static_libprimecount STATIC ${LIB_SRC})
    set_target_properties(static_libprimecount PROPERTIES OUTPUT_NAME primecount)
    target_link_libraries(static_libprimecount static_libprimesieve ${MPI_CXX_LIBRARIES})

    target_include_directories(static_libprimecount PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

    install(TARGETS static_libprimecount
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# primecount binary ##################################################

if(BUILD_PRIMECOUNT)
    add_executable(primecount ${BIN_SRC})
    target_link_libraries(primecount libprimecount)
    install(TARGETS primecount
            DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()

# Install header #####################################################

install(FILES include/primecount.hpp
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
