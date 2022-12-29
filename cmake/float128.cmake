# Check if the compiler and CPU architecture support
# 128-bit floating point numbers & libquadmath.

# As of 2022 only GCC supports libquadmath.

# Note that GCC currently warns when compiling using -pedantic
# because of the use of the non-standard Q suffix:
# "warning: non-standard suffix on floating constant [-Wpedantic]".
# And there is currently (2022) no way to suppress this warning,
# see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=92826.
# Because of this issue this test will fail when compiling with
# GCC and -pedantic -Werror (which is used e.g. in our CI tests).

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_LIBRARIES "quadmath")

# Code is similar to what's used in RiemannR.cpp
check_cxx_source_compiles("
    #include <quadmath.h>
    int main() {
        __float128 t = 10000000000.0Q;
        __float128 old_term = FLT128_MAX;

        for (int i = 0; i < 10; i++)
        {
            __float128 term = logq(t);
            if (fabsq(term) >= fabsq(old_term))
                break;
            t -= term;
            old_term = term;
        }

        return t > 0;
    }" float128)

cmake_pop_check_state()
