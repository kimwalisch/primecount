# Check for environment variable getter/setter functions.

if(WIN32)
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists(_putenv_s "stdlib.h" WIN_putenv_s)
    if(WIN_putenv_s)
        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("
            #include <stdlib.h>

            int main()
            {
                size_t required = 0;
                return getenv_s(&required, nullptr, 0, \"GETENV_TEST\");
            }" WIN_getenv_s)
    endif()
endif()

if(NOT WIN_getenv_s AND NOT WIN_putenv_s)
    check_cxx_source_compiles("
        #include <stdlib.h>

        int main()
        {
            return setenv(\"SETENV_TEST\", \"1\", 0);
        }" POSIX_setenv)
endif()
