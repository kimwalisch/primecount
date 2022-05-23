include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(OpenMP QUIET)

if(NOT (OpenMP_FOUND OR OpenMP_CXX_FOUND))
    message(STATUS "Performing Test OpenMP")
    message(STATUS "Performing Test OpenMP - Failed")
endif()

# Check if libatomic is required (usually by Clang)
# to support 128-bit integers with OpenMP.
if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

    if(TARGET OpenMP::OpenMP_CXX)
        set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")
    else()
        set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS}")
    endif()

    if(DISABLE_INT128)
        set(CMAKE_REQUIRED_DEFINITIONS "-D${DISABLE_INT128}")
    endif()

    # Check if compiles without libatomic
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
        find_library(LIB_ATOMIC NAMES atomic libatomic.so.1)

        if(LIB_ATOMIC)
            set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS}")
            set(CMAKE_REQUIRED_LIBRARIES "${LIB_ATOMIC}")

            # Check if compiles with libatomic
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
        endif()

        if(NOT OpenMP_with_libatomic)
            set(LIB_ATOMIC "")
        endif()
    endif()

    cmake_pop_check_state()

    # OpenMP has been tested successfully, enable it
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
