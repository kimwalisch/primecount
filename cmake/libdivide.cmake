include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

check_cxx_source_compiles("
    #include <libdivide.h>
    #include <stdint.h>
    int main() {
        libdivide::divider<uint64_t, libdivide::BRANCHFREE> d(3);
        uint64_t n = 1000000000;
        return (n / d == n / 3) ? 0 : 1;
    }" libdivide_branchfree)

cmake_pop_check_state()
