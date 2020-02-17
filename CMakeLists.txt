cmake_minimum_required(VERSION 3.4)
project(primecount CXX)
set(PRIMECOUNT_VERSION_MAJOR 5)
set(PRIMECOUNT_VERSION_MINOR 4)
set(PRIMECOUNT_VERSION "${PRIMECOUNT_VERSION_MAJOR}.${PRIMECOUNT_VERSION_MINOR}")

# Build options ######################################################

option(BUILD_PRIMECOUNT    "Build the primecount binary"           ON)
option(BUILD_LIBPRIMESIEVE "Build libprimesieve"                   ON)
option(BUILD_SHARED_LIBS   "Build the shared libprimecount"        OFF)
option(BUILD_STATIC_LIBS   "Build the static libprimecount"        ON)
option(BUILD_MANPAGE       "Regenerate man page using a2x program" OFF)
option(BUILD_TESTS         "Build the test programs"               OFF)

option(WITH_POPCNT         "Use the POPCNT instruction"            ON)
option(WITH_OPENMP         "Enable OpenMP multi-threading"         ON)
option(WITH_LIBDIVIDE      "Use libdivide.h"                       ON)
option(WITH_DIV32          "Use 32-bit division instead of 64-bit division whenever possible" ON)
option(WITH_FLOAT128       "Use __float128 (requires libquadmath)" OFF)
option(WITH_MPI            "Enable MPI support"                    OFF)

if(WIN32)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS)
    message(FATAL_ERROR "One or both of BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be set to ON")
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
            src/common.cpp
            src/FactorTable.cpp
            src/RiemannR.cpp
            src/P2.cpp
            src/P3.cpp
            src/PhiTiny.cpp
            src/PiTable.cpp
            src/SegmentedPiTable.cpp
            src/S1.cpp
            src/Sieve.cpp
            src/LoadBalancer.cpp
            src/S2Status.cpp
            src/generate.cpp
            src/nth_prime.cpp
            src/phi.cpp
            src/pi_legendre.cpp
            src/pi_lehmer.cpp
            src/pi_meissel.cpp
            src/pi_primesieve.cpp
            src/print.cpp
            src/lmo/pi_lmo1.cpp
            src/lmo/pi_lmo2.cpp
            src/lmo/pi_lmo3.cpp
            src/lmo/pi_lmo4.cpp
            src/lmo/pi_lmo5.cpp
            src/lmo/pi_lmo_parallel.cpp
            src/deleglise-rivat/S2_hard.cpp
            src/deleglise-rivat/S2_trivial.cpp
            src/deleglise-rivat/pi_deleglise_rivat3.cpp
            src/gourdon/pi_gourdon.cpp
            src/gourdon/Phi0.cpp
            src/gourdon/B.cpp
            src/gourdon/D4.cpp
            src/gourdon/Sigma.cpp)

# Required includes ##################################################

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)
include(CMakePushCheckState)

# Use 32-bit integer division ########################################

# Check at runtime if the dividend and divisor are < 2^32 and
# if so use 32-bit integer division instead of 64-bit integer
# division. On most CPUs before 2020 this significantly
# improves performance.
if(WITH_DIV32)
    set(ENABLE_DIV32 "ENABLE_DIV32")
endif()

# Check if libdivide.h compiles ######################################

if(WITH_LIBDIVIDE)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

    check_cxx_source_compiles("
        #include <libdivide.h>
        #include <stdint.h>
        int main() {
            libdivide::divider<uint64_t, libdivide::BRANCHFREE> d(3);
            uint64_t n = 1000000000;
            return (n / d == n / 3) ? 0 : 1;
        }" libdivide_branchfree)

    cmake_pop_check_state()
endif()

if(WITH_LIBDIVIDE AND libdivide_branchfree)
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

# Check for MPI (Message Passing Interface) ##########################

if(WITH_MPI)
    find_package(MPI REQUIRED)

    if(TARGET MPI::MPI_CXX)
        set(LIB_MPI "MPI::MPI_CXX")
    else()
        include_directories("${MPI_INCLUDE_PATH}")
    endif()

    set(ENABLE_MPI "ENABLE_MPI")
    set(LIB_SRC ${LIB_SRC}
                src/mpi/P2_mpi.cpp
                src/mpi/S2_hard_mpi.cpp
                src/mpi/MpiLoadBalancer.cpp
                src/mpi/MpiMsg.cpp)

    if(libdivide_branchfree)
        set(LIB_SRC ${LIB_SRC} src/mpi/S2_easy_mpi_libdivide.cpp)
    else()
        set(LIB_SRC ${LIB_SRC} src/mpi/S2_easy_mpi.cpp)
    endif()
endif()

# Silence GCC switch fall through warnings ###########################

check_cxx_compiler_flag(-Wno-implicit-fallthrough Wno_fallthrough)

if(Wno_fallthrough)
    set_source_files_properties(src/Sieve.cpp PROPERTIES COMPILE_FLAGS -Wno-implicit-fallthrough)
endif()

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    find_package(OpenMP QUIET)

    if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
        set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

        check_cxx_source_compiles("
            #include <int128_t.hpp>
            int main() {
                primecount::maxint_t sum = 0;
                #pragma omp parallel for reduction(+: sum)
                for (int i = 0; i < 100; i++)
                    sum += i;
                return 0;
            }" OpenMP)

        if(NOT OpenMP)
            find_library(LIB_ATOMIC NAMES atomic libatomic.so.1)

            if(NOT LIB_ATOMIC)
                set(LIB_ATOMIC "")
            else()
                set(CMAKE_REQUIRED_LIBRARIES "${LIB_ATOMIC}")

                check_cxx_source_compiles("
                    #include <int128_t.hpp>
                    int main() {
                        primecount::maxint_t sum = 0;
                        #pragma omp parallel for reduction(+: sum)
                        for (int i = 0; i < 100; i++)
                            sum += i;
                        return 0;
                    }" OpenMP_with_libatomic)

                if(NOT OpenMP_with_libatomic)
                    set(LIB_ATOMIC "")
                endif()
            endif()
        endif()

        cmake_pop_check_state()

        if(OpenMP OR OpenMP_with_libatomic)
            if(TARGET OpenMP::OpenMP_CXX)
                set(LIB_OPENMP "OpenMP::OpenMP_CXX")
            else()
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            endif()
        endif()
    else()
        message(STATUS "Performing Test OpenMP")
        message(STATUS "Performing Test OpenMP - Failed")
    endif()
endif()

# Find POPCNT compiler flag ##########################################

if(WITH_POPCNT)
    set(ENABLE_POPCNT "ENABLE_POPCNT")

    if("${CMAKE_CXX_COMPILER}" MATCHES "icpc")
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS -Werror)
        check_cxx_compiler_flag(-msse4.2 msse42)
        cmake_pop_check_state()

        if(msse42)
            set(POPCNT_FLAG "-msse4.2")
        endif()
    endif()

    if(NOT POPCNT_FLAG)
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS -Werror)
        check_cxx_compiler_flag(-mpopcnt mpopcnt)
        cmake_pop_check_state()

        if(mpopcnt)
            set(POPCNT_FLAG "-mpopcnt")
        endif()
    endif()
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

    set(BUILD_MANPAGE "${COPY_BUILD_MANPAGE}" CACHE BOOL "Regenerate man page using a2x" FORCE)
    set(BUILD_TESTS "${COPY_BUILD_TESTS}" CACHE BOOL "Build test programs" FORCE)
else()
    find_package(primesieve REQUIRED)
endif()

# Testing ############################################################

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# libprimecount ######################################################

if(BUILD_SHARED_LIBS)
    add_library(libprimecount SHARED ${LIB_SRC})
    set_target_properties(libprimecount PROPERTIES OUTPUT_NAME primecount)
    set_target_properties(libprimecount PROPERTIES SOVERSION ${PRIMECOUNT_VERSION_MAJOR})
    set_target_properties(libprimecount PROPERTIES VERSION ${PRIMECOUNT_VERSION})
    target_compile_definitions(libprimecount PRIVATE "${ENABLE_POPCNT}" "${ENABLE_DIV32}" "${ENABLE_MPI}")
    target_compile_options(libprimecount PRIVATE "${POPCNT_FLAG}")
    target_link_libraries(libprimecount PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_MPI}" "${LIB_QUADMATH}" "${LIB_ATOMIC}")

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
    target_compile_definitions(libprimecount-static PRIVATE "${ENABLE_POPCNT}" "${ENABLE_DIV32}" "${ENABLE_MPI}")
    target_compile_options(libprimecount-static PRIVATE "${POPCNT_FLAG}")
    target_link_libraries(libprimecount-static PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_MPI}" "${LIB_QUADMATH}" "${LIB_ATOMIC}")

    if(NOT BUILD_SHARED_LIBS)
        # make sure libprimecount is always defined
        add_library(libprimecount ALIAS libprimecount-static)
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

# primecount binary ##################################################

if(BUILD_PRIMECOUNT)
    add_executable(primecount ${BIN_SRC})
    target_compile_definitions(primecount PRIVATE "${ENABLE_MPI}")
    target_link_libraries(primecount PRIVATE libprimecount primesieve::primesieve "${LIB_MPI}")
    target_compile_features(primecount PRIVATE cxx_auto_type)
    install(TARGETS primecount DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Install header #####################################################

install(FILES include/primecount.h
              include/primecount.hpp
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
