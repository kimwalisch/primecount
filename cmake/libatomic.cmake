# This file first tries to find out whether it is necessary to
# add libatomic to the linker flags. Usually this is required
# on old 32-bit CPUs. Then it tries to find libatomic, ideally
# we would have liked to use an official CMake module to find
# libatomic, however there is none.
# See discussion: https://github.com/kimwalisch/primesieve/issues/141

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# <atomic> requires C++11 or later
if(NOT PRIMECOUNT_COMPILER_SUPPORTS_CPP11)
    if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    endif()
endif()

set(atomic64_TEST_SOURCE "
    #include <atomic>
    #include <stdint.h>
    int main() {
        std::atomic<uint64_t> x;
        x = 1;
        x.fetch_and(3);
        x--;
        return (int) x + (int) x.is_lock_free();
    }")

# Check if code compiles without libatomic.
# Should always work on CPUs >= 64-bits
check_cxx_source_compiles("${atomic64_TEST_SOURCE}" PRIMECOUNT_ATOMIC64)

# Our code requires libatomic to compile
if(NOT PRIMECOUNT_ATOMIC64)
    # First try -latomic as the compiler may find libraries
    # in directories which CMake does not search.
    set(CMAKE_REQUIRED_LIBRARIES "-latomic")
    check_cxx_source_compiles("${atomic64_TEST_SOURCE}" PRIMECOUNT_ATOMIC64_WITH_LATOMIC)

    if(PRIMECOUNT_ATOMIC64_WITH_LATOMIC)
        list(APPEND PRIMECOUNT_LINK_LIBRARIES "-latomic")
        string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "-latomic ")
    else()
        find_library(PRIMECOUNT_LIBATOMIC NAMES atomic atomic.so.1 libatomic.so.1)

        if(PRIMECOUNT_LIBATOMIC)
            set(CMAKE_REQUIRED_LIBRARIES "${PRIMECOUNT_LIBATOMIC}")
            check_cxx_source_compiles("${atomic64_TEST_SOURCE}" PRIMECOUNT_ATOMIC64_WITH_LIBATOMIC_PATH)
        endif()

        if(PRIMECOUNT_ATOMIC64_WITH_LIBATOMIC_PATH)
            list(APPEND PRIMECOUNT_LINK_LIBRARIES "${PRIMECOUNT_LIBATOMIC}")
        else()
            message(FATAL_ERROR "Failed to compile std::atomic, libatomic likely not found!")
        endif()
    endif()
endif()

cmake_pop_check_state()
