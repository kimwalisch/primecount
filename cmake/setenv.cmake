# Check for environment variable setter functions.

if(WIN32)
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists(_putenv_s "stdlib.h" WIN_putenv_s)
endif()

if(NOT WIN_putenv_s)
    include(CheckCXXSourceCompiles)

    check_cxx_source_compiles("
        #include <stdlib.h>

        int main()
        {
            return setenv(\"SETENV_TEST\", \"1\", 0);
        }" POSIX_setenv)
endif()
