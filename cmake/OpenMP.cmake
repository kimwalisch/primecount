# Check if OpenMP supports 128-bit integers out of the box
# or if we have to link against libatomic.

set(PRIMECOUNT_WITH_OPENMP OFF)

if(NOT WITH_OPENMP)
    return()
endif()

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(OpenMP QUIET)

if(NOT (TARGET OpenMP::OpenMP_CXX))
    message(STATUS "Performing Test PRIMECOUNT_OPENMP")
    message(STATUS "Performing Test PRIMECOUNT_OPENMP - Failed")
endif()

# CMake has found OpenMP, now we need to check
# if OpenMP supports 128-bit integers.
if(TARGET OpenMP::OpenMP_CXX)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")

    set(OpenMP_TEST_SOURCE "
        #include <int128_t.hpp>
        #include <omp.h>
        #include <stdint.h>
        #include <iostream>
        int main(int, char** argv) {
            using primecount::maxint_t;
            uintptr_t n = (uintptr_t) argv;
            maxint_t sum = (maxint_t) n;
            int iters = (int) n;
            #pragma omp parallel for reduction(+: sum)
            for (int i = 0; i < iters; i++)
                sum += (i / 3) * omp_get_thread_num();
            std::cout << (long) sum;
            return 0;
        }")

    # Our <int128_t.hpp> requires C++11 or later
    if(NOT PRIMECOUNT_COMPILER_SUPPORTS_CPP11)
        if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION}")
        elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
        endif()
    endif()

    if(DISABLE_INT128)
        set(CMAKE_REQUIRED_DEFINITIONS "-D${DISABLE_INT128}")
    endif()

    set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

    # Check if OpenMP supports 128-bit integers
    check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" PRIMECOUNT_OPENMP)

    if(NOT PRIMECOUNT_OPENMP)
        # Try if OpenMP works if we link against libatomic.
        # This is sometimes required for LLVM/Clang.
        # First try -latomic as the compiler may find libraries
        # in directories which CMake does not search.
        set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX" "-latomic")
        check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" PRIMECOUNT_OPENMP_WITH_LATOMIC)

        if(PRIMECOUNT_OPENMP_WITH_LATOMIC)
            list(APPEND PRIMECOUNT_LINK_LIBRARIES "-latomic")
            string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "-latomic ")
        else()
            find_library(PRIMECOUNT_LIB_ATOMIC NAMES atomic atomic.so.1 libatomic.so.1)

            if(PRIMECOUNT_LIB_ATOMIC)
                set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX" "${PRIMECOUNT_LIB_ATOMIC}")
                check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" PRIMECOUNT_OPENMP_WITH_LIBATOMIC_PATH)

                if(PRIMECOUNT_OPENMP_WITH_LIBATOMIC_PATH)
                    list(APPEND PRIMECOUNT_LINK_LIBRARIES "${PRIMECOUNT_LIB_ATOMIC}")
                endif()
            endif()
        endif()
    endif()

    cmake_pop_check_state()

    # OpenMP has been tested successfully, enable it
    if(PRIMECOUNT_OPENMP OR
       PRIMECOUNT_OPENMP_WITH_LATOMIC OR
       PRIMECOUNT_OPENMP_WITH_LIBATOMIC_PATH)
        set(PRIMECOUNT_WITH_OPENMP ON)
        list(APPEND PRIMECOUNT_LINK_LIBRARIES "OpenMP::OpenMP_CXX")

        # Create list of private libs for pkg-config/pkgconf
        foreach(X IN LISTS OpenMP_CXX_LIB_NAMES)
            string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "-l${X} ")
        endforeach()
    endif()
endif()

# If we are using LLVM OpenMP we check if the compiler
# supports setenv() to tune the LLVM OpenMP options.
if(PRIMECOUNT_WITH_OPENMP)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")
    check_cxx_source_compiles("
        #include <omp.h>
        int main() {
            kmp_set_blocktime(1);
            return 0;
        }" PRIMECOUNT_LLVM_OPENMP)
    cmake_pop_check_state()

    if(PRIMECOUNT_LLVM_OPENMP)
        include("${PROJECT_SOURCE_DIR}/cmake/setenv.cmake")
    endif()
endif()

# OpenMP is not supported, print warning message
if(NOT PRIMECOUNT_WITH_OPENMP)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|LLVM")
        message(WARNING "Install the OpenMP library (libomp) to enable multithreading in primecount!")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(WARNING "Install the OpenMP library (libgomp) to enable multithreading in primecount!")
    else()
        message(WARNING "Install the OpenMP library to enable multithreading in primecount!")
    endif()
endif()
