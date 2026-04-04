if(NOT DEFINED ASM_FILE)
    message(FATAL_ERROR "ASM_FILE is not set.")
endif()

string(REGEX REPLACE "^\"(.*)\"$" "\\1" ASM_FILE "${ASM_FILE}")

if(NOT EXISTS "${ASM_FILE}")
    message(FATAL_ERROR "Assembly file not found: ${ASM_FILE}")
endif()

file(READ "${ASM_FILE}" primecount_conditional_move_asm)
string(TOLOWER "${primecount_conditional_move_asm}" primecount_conditional_move_asm)

string(REGEX MATCHALL "(^|[\r\n])[ \t]*cmov[a-z0-9]+([ \t]|$)"
       primecount_cmov_matches
       "${primecount_conditional_move_asm}")
list(LENGTH primecount_cmov_matches primecount_cmov_count)

if(NOT primecount_cmov_count EQUAL 2)
    message(FATAL_ERROR
            "Expected exactly 2 cmov* instructions in ${ASM_FILE}, "
            "but found ${primecount_cmov_count}.")
endif()

message(STATUS "Found exactly 2 cmov* instructions in ${ASM_FILE}")
