cmake_minimum_required(VERSION 3.4...3.22)
project(primecount CXX)
set(PRIMECOUNT_VERSION_MAJOR 7)
set(PRIMECOUNT_VERSION_MINOR 6)
set(PRIMECOUNT_VERSION "${PRIMECOUNT_VERSION_MAJOR}.${PRIMECOUNT_VERSION_MINOR}")

# Build options ######################################################

option(BUILD_PRIMECOUNT    "Build the primecount binary"           ON)
option(BUILD_LIBPRIMESIEVE "Build libprimesieve"                   ON)
option(BUILD_SHARED_LIBS   "Build the shared libprimecount"        OFF)
option(BUILD_STATIC_LIBS   "Build the static libprimecount"        ON)
option(BUILD_MANPAGE       "Regenerate man page using a2x program" OFF)
option(BUILD_TESTS         "Build the test programs"               OFF)

option(WITH_POPCNT          "Use the POPCNT instruction"            ON)
option(WITH_LIBDIVIDE       "Use libdivide.h"                       ON)
option(WITH_OPENMP          "Enable OpenMP multi-threading"         ON)
option(WITH_MSVC_CRT_STATIC "Link primecount.lib with /MT instead of the default /MD" OFF)
option(WITH_FLOAT128        "Use __float128 (requires libquadmath), increases precision of Li(x) & RiemannR" OFF)
option(WITH_JEMALLOC        "Use jemalloc allocator"                OFF)

# When using WITH_DIV32=ON primecount checks at runtime
# if 32-bit division can be used instead of 64-bit
# division which improves performance on old x64 CPUs
# (<= 2020). However Clang performs this optimization by
# itself, so we can disable WITH_DIV32.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "LLVM")
    option(WITH_DIV32 "Use 32-bit division instead of 64-bit division whenever possible" OFF)
else()
    option(WITH_DIV32 "Use 32-bit division instead of 64-bit division whenever possible" ON)
endif()

# Option sanity checks ###############################################

# We don't yet support building libprimecount as a shared DLL
# library on Windows. If you need this let me know...
if(WIN32)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS)
    message(FATAL_ERROR "One or both of BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be set to ON!")
endif()

if(STATICALLY_LINK_LIBPRIMECOUNT AND NOT BUILD_STATIC_LIBS)
    message(FATAL_ERROR "Using STATICALLY_LINK_LIBPRIMECOUNT=ON requires BUILD_STATIC_LIBS=ON!")
endif()

# Static vs. shared linking ##########################################

# By default libprimecount is linked statically if only the static
# libprimecount has been built. If libprimecount has been built as
# both a static and shared library then the shared libprimecount will
# be used for linking.

# It is possible to force static linking by setting:
# cmake -DSTATICALLY_LINK_LIBPRIMECOUNT=ON
# This way libprimecount will be linked statically even if
# libprimecount has been built as both a static and shared library.

if(NOT STATICALLY_LINK_LIBPRIMECOUNT)
    if(BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
        set(STATICALLY_LINK_LIBPRIMECOUNT TRUE)
    endif()
endif()

# Set default build type to Release ##################################

if(NOT CMAKE_VERSION VERSION_LESS 3.9)
    get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
elseif(CMAKE_CONFIGURATION_TYPES)
    set(isMultiConfig TRUE)
endif()

if(NOT isMultiConfig AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# primecount binary source files #####################################

set(BIN_SRC src/app/cmdoptions.cpp
            src/app/main.cpp
            src/app/help.cpp
            src/app/test.cpp)

# primecount library source files ####################################

set(LIB_SRC src/api.cpp
            src/api_c.cpp
            src/BitSieve240.cpp
            src/FactorTable.cpp
            src/RiemannR.cpp
            src/P2.cpp
            src/P3.cpp
            src/PhiTiny.cpp
            src/PiTable.cpp
            src/S1.cpp
            src/Sieve.cpp
            src/LoadBalancerP2.cpp
            src/LoadBalancerS2.cpp
            src/StatusS2.cpp
            src/generate.cpp
            src/nth_prime.cpp
            src/phi.cpp
            src/pi_legendre.cpp
            src/pi_lehmer.cpp
            src/pi_meissel.cpp
            src/pi_primesieve.cpp
            src/print.cpp
            src/util.cpp
            src/lmo/pi_lmo1.cpp
            src/lmo/pi_lmo2.cpp
            src/lmo/pi_lmo3.cpp
            src/lmo/pi_lmo4.cpp
            src/lmo/pi_lmo5.cpp
            src/lmo/pi_lmo_parallel.cpp
            src/deleglise-rivat/S2_hard.cpp
            src/deleglise-rivat/S2_trivial.cpp
            src/deleglise-rivat/pi_deleglise_rivat.cpp
            src/gourdon/pi_gourdon.cpp
            src/gourdon/Phi0.cpp
            src/gourdon/B.cpp
            src/gourdon/D.cpp
            src/gourdon/LoadBalancerAC.cpp
            src/gourdon/SegmentedPiTable.cpp
            src/gourdon/Sigma.cpp)

# Required includes ##################################################

include(GNUInstallDirs)

# Use -Wno-uninitialized with GCC compiler ###########################

# GCC's -Wuninitialized enabled with -Wall -pedantic causes
# false postive warnings in libprimecount:
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=107287
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(WNO_UNINITIALIZED "-Wno-uninitialized")
endif()

# libprimecount compile definitions ##################################

# Check at runtime if the dividend and divisor are < 2^32 and
# if so use 32-bit integer division instead of 64-bit integer
# division. On most CPUs before 2020 this significantly
# improves performance.
if(WITH_DIV32)
    set(ENABLE_DIV32 "ENABLE_DIV32")
endif()

# Enable expensive debugging assertions.
# These assertions enable e.g. bounds checks for the pod_vector
# and pod_array types.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(ENABLE_ASSERT "ENABLE_ASSERT")
endif()

# Check if int128_t is supported #####################################

include("${PROJECT_SOURCE_DIR}/cmake/int128_t.cmake")

# Use libdivide.h (fast integer divison) #############################

if(WITH_LIBDIVIDE)
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
    set(LIB_SRC ${LIB_SRC} src/gourdon/AC_libdivide.cpp)
else()
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy.cpp)
    set(LIB_SRC ${LIB_SRC} src/gourdon/AC.cpp)
endif()

# Enable __float128 support (requires libquadmath) ###################

if(WITH_FLOAT128)
    set(LIB_QUADMATH "quadmath")
    set_source_files_properties(src/RiemannR.cpp PROPERTIES COMPILE_DEFINITIONS HAVE_FLOAT128)
endif()

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    include("${PROJECT_SOURCE_DIR}/cmake/OpenMP.cmake")
endif()

# Check if x86 CPU supports POPCNT instruction #######################

if(WITH_POPCNT)
    include("${PROJECT_SOURCE_DIR}/cmake/popcnt.cmake")
endif()

# libprimesieve ######################################################

# By default the libprimesieve dependency is built from source
# (BUILD_LIBPRIMESIEVE=ON). However when packaging primecount
# for e.g. a Linux distro this is not a good idea. For this use
# case it is better to install the libprimesieve package
# (e.g. libprimesieve-dev) and link against that version.

if(BUILD_LIBPRIMESIEVE)
    set(COPY_BUILD_MANPAGE "${BUILD_MANPAGE}")
    set(COPY_BUILD_TESTS "${BUILD_TESTS}")
    set(BUILD_MANPAGE OFF CACHE BOOL "Build primesieve manpage" FORCE)
    set(BUILD_TESTS OFF CACHE BOOL "Build primesieve tests" FORCE)
    option(BUILD_PRIMESIEVE "Build primesieve binary" OFF)

    add_subdirectory(lib/primesieve)

    if(BUILD_SHARED_LIBS)
        target_compile_options(libprimesieve PRIVATE "${POPCNT_FLAG}")
    endif()
    if(BUILD_STATIC_LIBS)
        target_compile_options(libprimesieve-static PRIVATE "${POPCNT_FLAG}")
    endif()

    set(BUILD_MANPAGE "${COPY_BUILD_MANPAGE}" CACHE BOOL "Regenerate man page using a2x" FORCE)
    set(BUILD_TESTS "${COPY_BUILD_TESTS}" CACHE BOOL "Build test programs" FORCE)
else()
    find_package(primesieve REQUIRED)

    if(primesieve_VERSION VERSION_LESS 11.0)
        message(FATAL_ERROR "Found libprimesive-${primesieve_VERSION}, but primecount requires libprimesive >= 11.0")
    endif()
endif()

# libprimecount ######################################################

if(BUILD_SHARED_LIBS)
    add_library(libprimecount SHARED ${LIB_SRC})
    set_target_properties(libprimecount PROPERTIES OUTPUT_NAME primecount)
    set_target_properties(libprimecount PROPERTIES SOVERSION ${PRIMECOUNT_VERSION_MAJOR})
    set_target_properties(libprimecount PROPERTIES VERSION ${PRIMECOUNT_VERSION})
    target_compile_options(libprimecount PRIVATE "${POPCNT_FLAG}" "${WNO_UNINITIALIZED}")
    target_compile_definitions(libprimecount PRIVATE "${DISABLE_INT128}" "${ENABLE_DIV32}" "${ENABLE_ASSERT}")
    target_link_libraries(libprimecount PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_QUADMATH}" "${LIB_ATOMIC}")

    target_compile_features(libprimecount
    PRIVATE
        cxx_constexpr
        cxx_lambdas
        cxx_range_for)

    target_include_directories(libprimecount PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

    install(TARGETS libprimecount
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# libprimecount-static ###############################################

if(BUILD_STATIC_LIBS)
    add_library(libprimecount-static STATIC ${LIB_SRC})
    set_target_properties(libprimecount-static PROPERTIES OUTPUT_NAME primecount)
    target_compile_options(libprimecount-static PRIVATE "${POPCNT_FLAG}" "${WNO_UNINITIALIZED}")
    target_compile_definitions(libprimecount-static PRIVATE "${DISABLE_INT128}" "${ENABLE_DIV32}" "${ENABLE_ASSERT}")
    target_link_libraries(libprimecount-static PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_QUADMATH}" "${LIB_ATOMIC}")

    if(WITH_MSVC_CRT_STATIC)
        set_target_properties(libprimecount-static PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()

    target_compile_features(libprimecount-static
    PRIVATE
        cxx_constexpr
        cxx_lambdas
        cxx_range_for)

    target_include_directories(libprimecount-static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

    install(TARGETS libprimecount-static
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Shared vs. static linking ##########################################

if(STATICALLY_LINK_LIBPRIMECOUNT)
    add_library(primecount::primecount ALIAS libprimecount-static)
else()
    add_library(primecount::primecount ALIAS libprimecount)
endif()

# primecount binary ##################################################

if(BUILD_PRIMECOUNT)
    add_executable(primecount ${BIN_SRC})
    target_link_libraries(primecount PRIVATE primecount::primecount primesieve::primesieve)
    target_compile_definitions(primecount PRIVATE "${DISABLE_INT128}" "${ENABLE_DIV32}" "${ENABLE_ASSERT}")
    target_compile_features(primecount PRIVATE cxx_auto_type)
    install(TARGETS primecount DESTINATION ${CMAKE_INSTALL_BINDIR})

    if(WITH_MSVC_CRT_STATIC)
        set_target_properties(primecount PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()
endif()

# Use jemalloc allocator #############################################

if(WITH_JEMALLOC)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JEMALLOC jemalloc)
    pkg_search_module(JEMALLOC REQUIRED jemalloc)

    if(BUILD_PRIMECOUNT)
        target_link_libraries(primecount PRIVATE ${JEMALLOC_LIBRARIES})
        target_include_directories(primecount PRIVATE ${JEMALLOC_INCLUDE_DIRS})
    endif()
    if(BUILD_SHARED_LIBS)
        target_link_libraries(libprimecount PRIVATE ${JEMALLOC_LIBRARIES})
        target_include_directories(libprimecount PRIVATE ${JEMALLOC_INCLUDE_DIRS})
    endif()
    if(BUILD_STATIC_LIBS)
        target_link_libraries(libprimecount-static PRIVATE ${JEMALLOC_LIBRARIES})
        target_include_directories(libprimecount-static PRIVATE ${JEMALLOC_INCLUDE_DIRS})
    endif()
endif()

# Install header #####################################################

install(FILES include/primecount.h
              include/primecount.hpp
        COMPONENT libprimecount-headers
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Regenerate man page ################################################

if(BUILD_MANPAGE)
    find_program(A2X_EXECUTABLE a2x)

    if(NOT A2X_EXECUTABLE)
        message(FATAL_ERROR "Missing a2x program (required for man page generation)!")
    else()
        message(STATUS "Found a2x: ${A2X_EXECUTABLE}")

        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${A2X_EXECUTABLE}
            ARGS --format=manpage
                 -D "${PROJECT_SOURCE_DIR}/doc"
                 "${PROJECT_SOURCE_DIR}/doc/primecount.txt"
            VERBATIM)
    endif()
endif()

# Install man page ###################################################

if(BUILD_PRIMECOUNT)
    install(FILES ${PROJECT_SOURCE_DIR}/doc/primecount.1
            DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()

# Install primecount.pc (pkg-config) #################################

configure_file(primecount.pc.in primecount.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/primecount.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Testing ############################################################

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
