include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# Our <FastDivider.hpp> requires C++11 or later
if(NOT PRIMECOUNT_COMPILER_SUPPORTS_CPP11)
    if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    endif()
endif()

# Use the same int128_t support as libprimecount
if(DISABLE_INT128)
    set(CMAKE_REQUIRED_DEFINITIONS "-D${DISABLE_INT128}")
endif()

set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

# FastDivider requires int128_t or MSVC's __umulh() intrinsic,
# hence it is not supported on 32-bit CPUs.
check_cxx_source_compiles("
    #include <FastDivider.hpp>
    #include <stdint.h>

    using namespace primecount;

    int main() {
        FastDivider divider(7);
        uint64_t n = 1000;
        return (n / divider == 142) ? 0 : 1;
    }" PRIMECOUNT_FASTDIVIDER)

if(PRIMECOUNT_FASTDIVIDER)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "ENABLE_FASTDIVIDER")
endif()

cmake_pop_check_state()
