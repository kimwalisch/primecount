cmake_minimum_required(VERSION 3.4...3.19)
project(primecount CXX)
set(PRIMECOUNT_VERSION_MAJOR 7)
set(PRIMECOUNT_VERSION_MINOR 2)
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
option(WITH_DIV32           "Use 32-bit division instead of 64-bit division whenever possible" ON)
option(WITH_MSVC_CRT_STATIC "Link primecount.lib with /MT instead of the default /MD" OFF)
option(WITH_FLOAT128        "Use __float128 (requires libquadmath)" OFF)
option(WITH_JEMALLOC        "Use jemalloc allocator"                OFF)

# libprimecount sanity checks ########################################

# We don't yet support building libprimecount as a shared DLL
# library on Windows. If you need this let me know...
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

        # We made sure OpenMP works. In 2021 OpenMP
        # does e.g. not work using Apple Clang.
        if(OpenMP OR OpenMP_with_libatomic)
            if(TARGET OpenMP::OpenMP_CXX)
                set(LIB_OPENMP "OpenMP::OpenMP_CXX")
            else()
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            endif()

            # Create list of OpenMP libs for pkg-config/pkgconf
            foreach(X IN LISTS OpenMP_CXX_LIB_NAMES)
                string(APPEND PKGCONFIG_LIBS_OPENMP "-l${X} ")
            endforeach()
        endif()
    else()
        message(STATUS "Performing Test OpenMP")
        message(STATUS "Performing Test OpenMP - Failed")
    endif()
endif()

# Find POPCNT compiler flag ##########################################

if(WITH_POPCNT)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS -Werror)
    check_cxx_compiler_flag(-mpopcnt mpopcnt)
    cmake_pop_check_state()

    if(mpopcnt)
        set(POPCNT_FLAG "-mpopcnt")
    endif()
else()
    # On the x86/x64 CPU architectures WITH_POPCNT=OFF disables the
    # POPCNT instruction and instead uses a portable integer popcount
    # algorithm. This deteriorates performance by up to 50%.
    # Disabling POPCNT on x86/x64 is necessary for old CPUs from
    # before 2008. On other CPU architectures POPCNT is always used
    # if it is available as this generally does not cause any issues.
    # Hence WITH_POPCNT=OFF is ignored on non x86/x64 CPUs.

    check_cxx_source_compiles("
        int main() {
            #if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
                return 0;
            #else
                Error: not MSVC compiler
            #endif
        }" MSVC_x86_or_x64)

    # We only need to disable POPCNT for MSVC on x86/x64.
    # For GCC/Clang __builtin_popcountll() will automatically
    # be portable if we compile without the -mpopcnt flag.
    if(MSVC_x86_or_x64)
        set(DISABLE_POPCNT "DISABLE_POPCNT")
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
    find_package(primesieve 7.6 REQUIRED)
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
    target_compile_definitions(libprimecount PRIVATE "${DISABLE_POPCNT}" "${ENABLE_DIV32}")
    target_compile_options(libprimecount PRIVATE "${POPCNT_FLAG}")
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
    target_compile_definitions(libprimecount-static PRIVATE "${DISABLE_POPCNT}" "${ENABLE_DIV32}")
    target_compile_options(libprimecount-static PRIVATE "${POPCNT_FLAG}")
    target_link_libraries(libprimecount-static PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_QUADMATH}" "${LIB_ATOMIC}")

    if(WITH_MSVC_CRT_STATIC)
        set_target_properties(libprimecount-static PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()

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
    target_link_libraries(primecount PRIVATE libprimecount primesieve::primesieve)
    target_compile_features(primecount PRIVATE cxx_auto_type)
    install(TARGETS primecount DESTINATION ${CMAKE_INSTALL_BINDIR})

    if(WITH_MSVC_CRT_STATIC)
        set_target_properties(primecount PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()
endif()

# Use jemalloc allocator (better scaling) ############################

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
