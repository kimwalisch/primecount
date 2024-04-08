include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# Check if we need to enable C++11 for our tests.
check_cxx_source_compiles("
    #include <limits>
    #include <type_traits>
    int main() {
        static_assert(std::numeric_limits<unsigned int>::min() == 0, \"\");
        static_assert(std::is_integral<int>::value, \"\");
        return 0;
    }" compiler_supports_cpp11)

if(NOT compiler_supports_cpp11)
    if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

# Test if the C++ compiler supports 128-bit integers.
# This test does not use the C++ standard library.
check_cxx_source_compiles("
    #include <int128_t.hpp>
    int main() {
        using primecount::int128_t;
        int128_t x = int128_t(1) << 100;
        x /= 123;
        return (x > 0) ? 0 : 1;
    }" int128)

if(int128)
    # Test if the C++ standard library supports
    # 128-bit integers out of the box.
    check_cxx_source_compiles("
        #include <int128_t.hpp>
        #include <algorithm>
        #include <limits>
        #include <type_traits>

        using namespace primecount;
        static_assert(std::numeric_limits<int128_t>::max() == int128_t((uint128_t(1) << 127) - 1), \"\");
        static_assert(std::numeric_limits<uint128_t>::max() == ~uint128_t(0), \"\");
        static_assert(std::numeric_limits<int128_t>::digits == 127, \"\");
        static_assert(std::numeric_limits<uint128_t>::digits == 128, \"\");
        static_assert(std::is_integral<int128_t>::value, \"\");
        static_assert(std::is_integral<uint128_t>::value, \"\");
        static_assert(std::is_signed<int128_t>::value, \"\");
        static_assert(std::is_unsigned<uint128_t>::value, \"\");
        static_assert(std::is_unsigned<typename std::make_unsigned<int128_t>::type>::value, \"\");

        int main() {
            int128_t x = 123;
            int128_t y = 1000;
            int128_t z = std::max(x, y);
            if (z == y)
                return 0;
            else
                return 1;
        }" int128_STL)

    if(NOT int128_STL)
        # Now let's try again if it works with our <int128_STL_patch.hpp>
        # header. This is required if the user compiles with
        # -std=c++11 instead of -std=gnu++11.
        check_cxx_source_compiles("
            #include <int128_t.hpp>
            #include <int128_STL_patch.hpp>
            #include <algorithm>
            #include <limits>
            #include <type_traits>

            using namespace primecount;
            static_assert(std::numeric_limits<int128_t>::max() == int128_t((uint128_t(1) << 127) - 1), \"\");
            static_assert(std::numeric_limits<uint128_t>::max() == ~uint128_t(0), \"\");
            static_assert(std::numeric_limits<int128_t>::digits == 127, \"\");
            static_assert(std::numeric_limits<uint128_t>::digits == 128, \"\");
            static_assert(std::is_integral<int128_t>::value, \"\");
            static_assert(std::is_integral<uint128_t>::value, \"\");
            static_assert(std::is_signed<int128_t>::value, \"\");
            static_assert(std::is_unsigned<uint128_t>::value, \"\");
            static_assert(std::is_unsigned<typename std::make_unsigned<int128_t>::type>::value, \"\");

            int main() {
                int128_t x = 123;
                int128_t y = 1000;
                int128_t z = std::max(x, y);
                if (z == y)
                    return 0;
                else
                    return 1;
            }" int128_STL_patch)

        if(int128_STL_patch)
            set(ENABLE_INT128_STL_PATCH "ENABLE_INT128_STL_PATCH")
        endif()
    endif()
endif()

if(int128 AND NOT int128_STL AND NOT int128_STL_patch)
    set(DISABLE_INT128 "DISABLE_INT128")
    message(WARNING "Your C++ standard library does not support 128-bit integers, "
                    "hence primecount will only support numbers <= 2^63.")
endif()

cmake_pop_check_state()
