include(CheckCXXSourceRuns)
include(CMakePushCheckState)
include(CheckCXXCompilerFlag)

cmake_push_check_state()
set(CMAKE_REQUIRED_FLAGS -Werror)
check_cxx_compiler_flag(-mpopcnt mpopcnt)
cmake_pop_check_state()

if(mpopcnt)
    set(POPCNT_FLAG "-mpopcnt")
endif()

# Since the POPCNT instruction is very important for performance
# (up to 50% speedup), we will only disable POPCNT if we can
# prove (using 2 tests) that the CPU does not support it.
cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")
set(CMAKE_REQUIRED_QUIET TRUE)

if(NOT CMAKE_CROSSCOMPILING)
    # GCC/Clang & x86 CPU
    if(mpopcnt)
        check_cxx_source_runs("
            #include <popcnt.hpp>
            #include <stdint.h>
            #include <iostream>
            int main(int, char** argv) {
                uintptr_t n = (uintptr_t) argv;
                std::cout << popcnt64((uint64_t) n);
                return 0;
            }" popcnt64_without_mpopcnt)

        if(popcnt64_without_mpopcnt)
            set(CMAKE_REQUIRED_FLAGS "${POPCNT_FLAG}")
            set(CMAKE_REQUIRED_QUIET FALSE)

            check_cxx_source_runs("
                #include <popcnt.hpp>
                #include <stdint.h>
                #include <iostream>
                int main(int, char** argv) {
                    uintptr_t n = (uintptr_t) argv;
                    std::cout << popcnt64((uint64_t) n);
                    return 0;
                }" cpu_supports_popcnt)

            if(NOT cpu_supports_popcnt)
                set(POPCNT_FLAG "")
            endif()
        endif()
    endif()
endif()

cmake_pop_check_state()
