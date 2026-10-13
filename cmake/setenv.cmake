# Check for environment variable getter/setter functions.

if(WIN32)
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists(_putenv_s "stdlib.h" PRIMECOUNT_WIN_PUTENV_S)
    if(PRIMECOUNT_WIN_PUTENV_S)
        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("
            #include <stdlib.h>

            int main()
            {
                size_t required = 0;
                return getenv_s(&required, nullptr, 0, \"GETENV_TEST\");
            }" PRIMECOUNT_WIN_GETENV_S)
    endif()
endif()

if(NOT PRIMECOUNT_WIN_GETENV_S AND NOT PRIMECOUNT_WIN_PUTENV_S)
    check_cxx_source_compiles("
        #include <stdlib.h>

        int main()
        {
            return setenv(\"SETENV_TEST\", \"1\", 0);
        }" PRIMECOUNT_POSIX_SETENV)
endif()
