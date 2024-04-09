include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# Our <int128_t.hpp> requires C++11 or later
if(NOT compiler_supports_cpp11)
    if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

# Our code requires very little 128-bit support from the C++ standard
# library because our <int128_t.hpp> contains portable replacements
# for the <type_traits> and <limits> headers. Our code basically only
# requires std::min() and std::max() support for int128_t.
check_cxx_source_compiles("
    #include <int128_t.hpp>
    #include <algorithm>

    using namespace primecount;

    int main() {
        int128_t x = int128_t(1) << 100;
        int128_t y = 1000;
        x /= 123;

        if (std::min(x, y) != y)
            return 1;

        long double z = 1001.59867;
        int128_t iz = (int128_t) z;

        if (std::max(y, iz) != iz)
            return 1;

        return 0;
    }" int128)

if(NOT int128)
    set(DISABLE_INT128 "DISABLE_INT128")
endif()

cmake_pop_check_state()
