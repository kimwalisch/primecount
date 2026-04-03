if(NOT DEFINED ASM_FILE)
    message(FATAL_ERROR "ASM_FILE is not set.")
endif()

string(REGEX REPLACE "^\"(.*)\"$" "\\1" ASM_FILE "${ASM_FILE}")

if(NOT EXISTS "${ASM_FILE}")
    message(FATAL_ERROR "Assembly file not found: ${ASM_FILE}")
endif()

file(READ "${ASM_FILE}" primecount_conditional_move_asm)
string(TOLOWER "${primecount_conditional_move_asm}" primecount_conditional_move_asm)

if(NOT primecount_conditional_move_asm MATCHES "(^|[\r\n])[ \t]*cmov[a-z0-9]+([ \t]|$)")
    message(FATAL_ERROR
            "Expected at least one cmov* instruction in ${ASM_FILE}, but none was found.")
endif()

message(STATUS "Found a cmov* instruction in ${ASM_FILE}")
