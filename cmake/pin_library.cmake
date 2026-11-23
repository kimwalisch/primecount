# Checks whether libprimecount can pin itself, i.e. prevent being
# unloaded, e.g. by dlclose(). The thread pool's workers run the
# library's code until the process exits, see ThreadPool.cpp.
# dladdr(), RTLD_NOLOAD and RTLD_NODELETE are not part of POSIX,
# but they are supported by glibc, musl, macOS and the BSDs. On
# Windows, ThreadPool.cpp uses GetModuleHandleExW() instead, which
# requires no check.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# Same code as pin_library() in ThreadPool.cpp, hence
# ThreadPool.cpp compiles whenever this check succeeds.
set(pin_library_TEST_SOURCE "
    #include <dlfcn.h>
    int main() {
        static char address;
        Dl_info info;
        if (dladdr(&address, &info) && info.dli_fname)
            dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE);
        return 0;
    }")

# glibc >= 2.34, musl, macOS and the BSDs
# provide these functions in libc.
check_cxx_source_compiles("${pin_library_TEST_SOURCE}" PRIMECOUNT_PIN_LIBRARY)

# Older glibc versions require libdl
if(NOT PRIMECOUNT_PIN_LIBRARY AND CMAKE_DL_LIBS)
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_DL_LIBS}")
    check_cxx_source_compiles("${pin_library_TEST_SOURCE}" PRIMECOUNT_PIN_LIBRARY_WITH_DL_LIBS)

    if(PRIMECOUNT_PIN_LIBRARY_WITH_DL_LIBS)
        list(APPEND PRIMECOUNT_LINK_LIBRARIES ${CMAKE_DL_LIBS})

        # CMAKE_DL_LIBS is a list of library names or linker
        # flags, e.g. "dl" on Linux and "-lld" on AIX.
        foreach(lib IN LISTS CMAKE_DL_LIBS)
            if(lib MATCHES "^-")
                string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "${lib} ")
            else()
                string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "-l${lib} ")
            endif()
        endforeach()
    endif()
endif()

if(PRIMECOUNT_PIN_LIBRARY OR PRIMECOUNT_PIN_LIBRARY_WITH_DL_LIBS)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "HAVE_PIN_LIBRARY")
endif()

cmake_pop_check_state()
