# Try if double floating-point division can be used instead of integer
# division without loss of precision if the dividend & divisor <= 2^53.
# Floating-point division if often much faster than integer division.
# This should work on virtually all CPU architectures, since
# virtually all CPU architectures use the IEEE 754 standard for
# floating-point arithmetic.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

if(CMAKE_CXX14_EXTENSION_COMPILE_OPTION)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX14_EXTENSION_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_CXX14_STANDARD_COMPILE_OPTION)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX14_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_FLAGS}")
endif()

check_cxx_source_compiles("
    #define ENABLE_DOUBLE_INTEGER_DIVISION
    #include <fast_div.hpp>

    int main() {
        // If the double mantissa has 53 bits we assume double floating-point
        // division is faster than 64-bit integer division.
        static_assert(std::numeric_limits<double>::digits == 53, \"double type mantissa bits != 53\");

        static_assert(fast_div(uint64_t(1ull << 53), uint32_t(1)) == 1ull << 53,
            \"fast_div(2^53, 1) failed using double floating-point division!\");
        static_assert(fast_div(uint64_t(1ull << 53), uint32_t(2)) == 1ull << 52,
            \"fast_div(2^53, 2) failed using double floating-point division!\");

        static_assert(fast_div(uint64_t(1ull << 53), uint64_t(1ull << 53)) == 1ull,
            \"fast_div(2^53, 2^53) failed using double floating-point division!\");
        static_assert(fast_div(uint64_t(1ull << 52), uint64_t(1ull << 53)) == 0ull,
            \"fast_div(2^52, 2^53) failed using double floating-point division!\");

        static_assert(fast_div(uint64_t((1ull << 53) - 1), uint32_t(2)) == 4503599627370495ull,
            \"fast_div(2^53-1, 2) failed using double floating-point division!\");
        static_assert(fast_div(uint64_t((1ull << 53) - 3), uint32_t(7)) == 1286742750677284ull,
            \"fast_div(2^53-3, 7) failed using double floating-point division!\");

        static_assert(fast_div(uint64_t(1ull << 53), uint64_t((1ull << 40) - 1)) == 8192ull,
            \"fast_div(2^53, 2^40-1) failed using double floating-point division!\");
        static_assert(fast_div(uint64_t((1ull << 53) - 1), uint64_t((1ull << 40) - 1)) == 8192ull,
            \"fast_div(2^53-1, 2^40-1) failed using double floating-point division!\");
        static_assert(fast_div(uint64_t((1ull << 53) - 10000), uint64_t((1ull << 40) - 1)) == 8191ull,
            \"fast_div(2^53-10000, 2^40-1) failed using double floating-point division!\");

        return 0;
    }" double_integer_division)

cmake_pop_check_state()
