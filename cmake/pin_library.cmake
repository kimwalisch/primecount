# Checks whether libprimecount can pin itself, i.e. prevent being
# unloaded, e.g. by dlclose(). The thread pool's workers run the
# library's code until the process exits, see ThreadPool.cpp.
# dladdr(), RTLD_NOLOAD and RTLD_NODELETE are not part of POSIX,
# but they are supported by glibc, musl, macOS and the BSDs.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

set(pin_library_TEST_SOURCE "
    #include <dlfcn.h>
    int main() {
        static char address;
        Dl_info info;
        if (dladdr(&address, &info) && info.dli_fname)
            return dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE) ? 0 : 1;
        return 1;
    }")

# glibc >= 2.34, musl, macOS and the BSDs
# provide these functions in libc.
check_cxx_source_compiles("${pin_library_TEST_SOURCE}" PRIMECOUNT_PIN_LIBRARY)

# Older glibc versions require libdl
if(NOT PRIMECOUNT_PIN_LIBRARY AND CMAKE_DL_LIBS)
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_DL_LIBS}")
    check_cxx_source_compiles("${pin_library_TEST_SOURCE}" PRIMECOUNT_PIN_LIBRARY_WITH_LIBDL)

    if(PRIMECOUNT_PIN_LIBRARY_WITH_LIBDL)
        list(APPEND PRIMECOUNT_LINK_LIBRARIES "${CMAKE_DL_LIBS}")
        string(APPEND PRIMECOUNT_PKGCONFIG_LIBS_PRIVATE "-l${CMAKE_DL_LIBS} ")
    endif()
endif()

if(PRIMECOUNT_PIN_LIBRARY OR PRIMECOUNT_PIN_LIBRARY_WITH_LIBDL)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "HAVE_PIN_LIBRARY")
endif()

cmake_pop_check_state()
