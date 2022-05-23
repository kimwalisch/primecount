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
    # MSVC compiler & x86 CPU
    else()
        set(CMAKE_REQUIRED_DEFINITIONS "-DDISABLE_POPCNT")
        check_cxx_source_runs("
            #include <popcnt.hpp>
            #include <stdint.h>
            #include <iostream>
            int main(int, char** argv) {
                #if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
                    uintptr_t n = (uintptr_t) argv;
                    std::cout << popcnt64((uint64_t) n);
                    return 0;
                #else
                    Error: not MSVC compiler
                #endif
            }" MSVC_without_popcnt)

        if(MSVC_without_popcnt)
            set(CMAKE_REQUIRED_DEFINITIONS "")
            set(CMAKE_REQUIRED_QUIET FALSE)

            check_cxx_source_runs("
                #include <popcnt.hpp>
                #include <stdint.h>
                #include <iostream>
                int main(int, char** argv) {
                    #if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
                        uintptr_t n = (uintptr_t) argv;
                        std::cout << popcnt64((uint64_t) n);
                        return 0;
                    #else
                        Error: not MSVC compiler
                    #endif
                }" cpu_supports_popcnt)

            if(NOT cpu_supports_popcnt)
                set(DISABLE_POPCNT "DISABLE_POPCNT")
            endif()
        endif()
    endif()
endif()

cmake_pop_check_state()
