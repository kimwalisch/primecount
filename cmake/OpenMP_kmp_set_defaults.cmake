# Check if the OpenMP runtime supports kmp_set_defaults().

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

if(TARGET OpenMP::OpenMP_CXX)
    set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")
endif()

check_cxx_source_compiles("
    #include <omp.h>

    int main()
    {
        kmp_set_defaults(\"OMP_WAIT_POLICY=ACTIVE\");
        return 0;
    }" OpenMP_kmp_set_defaults)

cmake_pop_check_state()
