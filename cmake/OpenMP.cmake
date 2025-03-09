# Check if OpenMP supports 128-bit integers out of the box
# or if we have to link against libatomic.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(OpenMP QUIET)

if(NOT (OpenMP_FOUND OR OpenMP_CXX_FOUND))
    message(STATUS "Performing Test OpenMP")
    message(STATUS "Performing Test OpenMP - Failed")
endif()

# CMake has found OpenMP, now we need to check
# if OpenMP supports 128-bit integers.
if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
    cmake_push_check_state()

    if(TARGET OpenMP::OpenMP_CXX)
        set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")
    else()
        set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS}")
    endif()

    # Our <int128_t.hpp> requires C++11 or later
    if(NOT compiler_supports_cpp11)
        if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${CMAKE_CXX11_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
        elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
        endif()
    endif()

    if(DISABLE_INT128)
        set(CMAKE_REQUIRED_DEFINITIONS "-D${DISABLE_INT128}")
    endif()

    set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

    # Check if OpenMP supports 128-bit integers
    check_cxx_source_compiles("
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
        }" OpenMP)

    if(NOT OpenMP)
        # Try if OpenMP works if we link against libatomic.
        # This is sometimes required for LLVM/Clang.
        find_library(LIB_ATOMIC NAMES atomic atomic.so.1 libatomic.so.1)

        if(NOT LIB_ATOMIC)
            # Some package managers like homebrew and macports store the compiler's
            # libraries in a subdirectory of the library directory. E.g. GCC
            # installed via homebrew stores libatomic at lib/gcc/13/libatomic.dylib
            # instead of lib/libatomic.dylib. CMake's find_library() cannot easily
            # be used to recursively find libraries. Therefore we use this workaround
            # here (try adding -latomic to linker options) for this use case.
            set(LIB_ATOMIC "-latomic")
        endif()

        set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}" "${LIB_ATOMIC}")

        # Check if OpenMP compiles with libatomic
        check_cxx_source_compiles("
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
            }" OpenMP_with_libatomic)

        if(OpenMP_with_libatomic)
            list(APPEND PRIMECOUNT_LINK_LIBRARIES "${LIB_ATOMIC}")
        else()
            set(LIB_ATOMIC "")
        endif()
    endif()

    cmake_pop_check_state()

    # OpenMP has been tested successfully, enable it
    if(OpenMP OR OpenMP_with_libatomic)
        if(TARGET OpenMP::OpenMP_CXX)
            list(APPEND PRIMECOUNT_LINK_LIBRARIES "OpenMP::OpenMP_CXX")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()

        # Create list of OpenMP libs for pkg-config/pkgconf
        foreach(X IN LISTS OpenMP_CXX_LIB_NAMES)
            string(APPEND PKGCONFIG_LIBS_OPENMP "-l${X} ")
        endforeach()
    endif()
endif()

# OpenMP test has failed, print warning message
if(NOT OpenMP AND NOT OpenMP_with_libatomic)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|LLVM")
        message(WARNING "Install the OpenMP library (libomp) to enable multithreading in primecount!")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(WARNING "Install the OpenMP library (libgomp) to enable multithreading in primecount!")
    else()
        message(WARNING "Install the OpenMP library to enable multithreading in primecount!")
    endif()
endif()
