# On x86 CPUs we need to enable the use of cpuid.cpp.
# If cpuid.cpp compiles we assume it is a x86 CPU.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}")

check_cxx_source_compiles("
    #include <src/arch/x86/cpuid.cpp>
    #include <iostream>

    int main()
    {
        int abcd[4];
        primecount::run_cpuid(1, 0, abcd);
        int ecx = abcd[2];

        if (ecx & (1 << 23)) == (1 << 23))
            std::cout << \"CPU supports POPCNT!\" << std::endl;
        else
            std::cout << \"CPU does not support POPCNT!\" << std::endl;

        return 0;
    }
" x86_cpuid)

cmake_pop_check_state()
