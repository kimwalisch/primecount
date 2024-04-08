include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(OpenMP QUIET)

if(NOT (OpenMP_FOUND OR OpenMP_CXX_FOUND))
    message(STATUS "Performing Test OpenMP")
    message(STATUS "Performing Test OpenMP - Failed")
endif()

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
        if (NOT DISABLE_INT128)
            # Check if OpenMP supports int128_t if we include
            # our <int128_OpenMP_patch.hpp> header.
            check_cxx_source_compiles("
                #include <int128_t.hpp>
                #include <int128_OpenMP_patch.hpp>
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
                }" OpenMP_int128_patch)
        endif()

        if(NOT OpenMP_int128_patch)
            # Finally try if OpenMP works if we link against libatomic.
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

            if(NOT OpenMP_with_libatomic)
                set(LIB_ATOMIC "")
            endif()
        endif()
    endif()

    cmake_pop_check_state()

    # OpenMP has been tested successfully, enable it
    if(OpenMP OR OpenMP_int128_patch OR OpenMP_with_libatomic)
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
if(NOT OpenMP AND NOT OpenMP_int128_patch AND NOT OpenMP_with_libatomic)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|LLVM")
        message(WARNING "Install the OpenMP library (libomp) to enable multithreading in primecount!")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(WARNING "Install the OpenMP library (libgomp) to enable multithreading in primecount!")
    else()
        message(WARNING "Install the OpenMP library to enable multithreading in primecount!")
    endif()
endif()
