cmake_minimum_required(VERSION 2.8.4)
project(primecount)
set(LIBPRIMECOUNT_VERSION "1.0.0")
set(LIBPRIMECOUNT_SOVERSION 1)
include_directories(include primesieve-latest/include)
set(CMAKE_BUILD_TYPE Release)

# Build options ######################################################

option(ENABLE_POPCNT "Enable POPCNT instruction" ON)
option(ENABLE_MPI "Enable MPI (Message Passing Interface)" OFF)

# Include Check* #####################################################

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)
include(CheckTypeSize)

# primecount binary source files #####################################

set(BIN_SRC src/app/cmdoptions.cpp
            src/app/help.cpp
            src/app/main.cpp)

# primecount library source files ####################################

set(LIB_SRC src/BitSieve.cpp
            src/FactorTable.cpp
            src/Li.cpp
            src/P2.cpp
            src/P3.cpp
            src/PhiTiny.cpp
            src/PiTable.cpp
            src/S1.cpp
            src/S2LoadBalancer.cpp
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
            src/lmo/pi_lmo_parallel1.cpp
            src/lmo/pi_lmo_parallel2.cpp
            src/lmo/pi_lmo_parallel3.cpp
            src/deleglise-rivat/S2_hard.cpp
            src/deleglise-rivat/S2_trivial.cpp
            src/deleglise-rivat/pi_deleglise_rivat1.cpp
            src/deleglise-rivat/pi_deleglise_rivat2.cpp
            src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp
            src/deleglise-rivat/pi_deleglise_rivat_parallel2.cpp)

# Check for MPI (Message Passing Interface) ##########################

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
    add_definitions(-DHAVE_MPI)

    set(LIB_SRC ${LIB_SRC} src/mpi/P2_mpi.cpp
                           src/mpi/deleglise-rivat
                           src/mpi/phi_vector.cpp
                           src/mpi/deleglise-rivat/S2_hard_mpi.cpp
                           src/mpi/deleglise-rivat/S2_hard_mpi_msg.cpp
                           src/mpi/deleglise-rivat/S2_hard_mpi_LoadBalancer.cpp)
endif()

# Check for OpenMP ###################################################

find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Check for -mpopcnt compiler option #################################

if(ENABLE_POPCNT)
    CHECK_CXX_COMPILER_FLAG(-mpopcnt MPOPCNT_FLAG)
    if(MPOPCNT_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif()
endif()

# Check for __int128_t ###############################################

check_type_size("__int128_t" INT128_T)
if(${INT128_T})
    add_definitions(-DHAVE___INT128_T)
    list(APPEND LIB_SRC src/deleglise-rivat/pi_deleglise_rivat_parallel3.cpp)
endif()

# Check if libdivide works ###########################################

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_SOURCE_DIR}/include)

check_cxx_source_runs("
    #include <libdivide.h>
    #include <stdint.h>
    #include <limits>

    int main()
    {
        uint64_t tests[20] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            63, 101, 511,
            1 << 5, 1 << 9, 1 << 20,
            std::numeric_limits<int32_t>::max(),
            std::numeric_limits<uint32_t>::max(),
            std::numeric_limits<int64_t>::max(),
            std::numeric_limits<uint64_t>::max() };

        for (int i = 0; i < 20; i++)
        {
            uint64_t x = tests[i];

            for (int j = 2; j < 1000; j++)
            {
                libdivide::divider<uint64_t, libdivide::BRANCHFREE> divider(j);
                if (x / j != x / divider)
                    return 1;
            }
        }

        return 0;
    }"
    HAVE_LIBDIVIDE)

# Include S2_easy* source files ######################################

if(HAVE_LIBDIVIDE)
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
    if(ENABLE_MPI)
        set(LIB_SRC ${LIB_SRC} src/mpi/deleglise-rivat/S2_easy_mpi_libdivide.cpp)
    endif()
else()
    set(LIB_SRC ${LIB_SRC} src/deleglise-rivat/S2_easy.cpp)
    if(ENABLE_MPI)
        set(LIB_SRC ${LIB_SRC} src/mpi/deleglise-rivat/S2_easy_mpi.cpp)
    endif()
endif()

# Find primesieve library ############################################

find_library(LIBPRIMESIEVE NAMES libprimesieve.a primesieve
                           PATHS ${CMAKE_SOURCE_DIR}/primesieve-latest/.libs)

# Build static library ###############################################

add_library(static_libprimecount STATIC ${LIB_SRC})
set_target_properties(static_libprimecount PROPERTIES OUTPUT_NAME primecount)

# Build shared library ###############################################

if(NOT (WIN32 OR MINGW OR CYGWIN))
    add_library(shared_libprimecount SHARED ${LIB_SRC})
    set_target_properties(shared_libprimecount PROPERTIES OUTPUT_NAME primecount)
    set_target_properties(shared_libprimecount PROPERTIES SOVERSION ${LIBPRIMECOUNT_SOVERSION})
    set_target_properties(shared_libprimecount PROPERTIES VERSION ${LIBPRIMECOUNT_VERSION})
    target_link_libraries(shared_libprimecount ${LIBPRIMESIEVE})
    install(TARGETS shared_libprimecount DESTINATION lib)
endif()

# Build binary #######################################################

add_executable(primecount ${BIN_SRC})
target_link_libraries(primecount static_libprimecount
                                 ${LIBPRIMESIEVE}
                                 ${MPI_CXX_LIBRARIES})

# Install primecount #################################################

install(TARGETS primecount DESTINATION bin)
install(TARGETS static_libprimecount DESTINATION lib)
install(FILES include/primecount.hpp DESTINATION include)
