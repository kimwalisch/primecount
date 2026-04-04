if(DEFINED PRIMECOUNT_CONFIGURE_CODEGEN_TEST)
    if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
        return()
    endif()

    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" PRIMECOUNT_SYSTEM_PROCESSOR_LOWER)
    string(TOLOWER "${CMAKE_CXX_COMPILER_ARCHITECTURE_ID}" PRIMECOUNT_COMPILER_ARCH_LOWER)

    if(NOT PRIMECOUNT_SYSTEM_PROCESSOR_LOWER MATCHES "^(x86_64|amd64)$" AND
       NOT PRIMECOUNT_COMPILER_ARCH_LOWER MATCHES "^(x64|x86_64|amd64)$")
        return()
    endif()

    get_filename_component(PRIMECOUNT_CODEGEN_NAME "${CMAKE_CURRENT_LIST_FILE}" NAME_WE)
    set(PRIMECOUNT_CODEGEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
    set(PRIMECOUNT_CODEGEN_SRC "${CMAKE_CURRENT_LIST_DIR}/${PRIMECOUNT_CODEGEN_NAME}.cpp")
    set(PRIMECOUNT_CODEGEN_CHECKER "${CMAKE_CURRENT_LIST_FILE}")
    set(PRIMECOUNT_CODEGEN_HEADERS
        "${PROJECT_SOURCE_DIR}/include/macros.hpp"
        "${PROJECT_SOURCE_DIR}/include/popcnt.hpp")
    set(PRIMECOUNT_CODEGEN_OBJ
        "${PRIMECOUNT_CODEGEN_DIR}/${PRIMECOUNT_CODEGEN_NAME}.obj")
    set(PRIMECOUNT_CODEGEN_ASM
        "${PRIMECOUNT_CODEGEN_DIR}/${PRIMECOUNT_CODEGEN_NAME}.s")
    set(PRIMECOUNT_CODEGEN_FLAGS "${CMAKE_CXX_FLAGS}")
    separate_arguments(PRIMECOUNT_CODEGEN_FLAGS NATIVE_COMMAND "${PRIMECOUNT_CODEGEN_FLAGS}")

    if(MSVC)
        set(PRIMECOUNT_CODEGEN_ASM
            "${PRIMECOUNT_CODEGEN_DIR}/${PRIMECOUNT_CODEGEN_NAME}.asm")
        add_custom_command(
            OUTPUT "${PRIMECOUNT_CODEGEN_ASM}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${PRIMECOUNT_CODEGEN_DIR}"
            COMMAND "${CMAKE_CXX_COMPILER}"
                    /nologo
                    /arch:AVX
                    /O2
                    /FAs
                    /c
                    ${PRIMECOUNT_CODEGEN_FLAGS}
                    "/I${PROJECT_SOURCE_DIR}/include"
                    "/Fo${PRIMECOUNT_CODEGEN_OBJ}"
                    "/Fa${PRIMECOUNT_CODEGEN_ASM}"
                    "${PRIMECOUNT_CODEGEN_SRC}"
            DEPENDS "${PRIMECOUNT_CODEGEN_SRC}" ${PRIMECOUNT_CODEGEN_HEADERS}
            COMMENT "Generating assembly for ${PRIMECOUNT_CODEGEN_NAME}"
            VERBATIM)
    else()
        add_custom_command(
            OUTPUT "${PRIMECOUNT_CODEGEN_ASM}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${PRIMECOUNT_CODEGEN_DIR}"
            COMMAND "${CMAKE_CXX_COMPILER}"
                    -mpopcnt
                    -O2
                    -S
                    -c
                    ${PRIMECOUNT_CODEGEN_FLAGS}
                    "-I${PROJECT_SOURCE_DIR}/include"
                    -o
                    "${PRIMECOUNT_CODEGEN_ASM}"
                    "${PRIMECOUNT_CODEGEN_SRC}"
            DEPENDS "${PRIMECOUNT_CODEGEN_SRC}" ${PRIMECOUNT_CODEGEN_HEADERS}
            COMMENT "Generating assembly for ${PRIMECOUNT_CODEGEN_NAME}"
            VERBATIM)
    endif()

    add_custom_target(${PRIMECOUNT_CODEGEN_NAME} ALL
        DEPENDS "${PRIMECOUNT_CODEGEN_ASM}")

    add_test(
        NAME ${PRIMECOUNT_CODEGEN_NAME}
        COMMAND ${CMAKE_COMMAND}
                "-DASM_FILE=${PRIMECOUNT_CODEGEN_ASM}"
                -P "${PRIMECOUNT_CODEGEN_CHECKER}")

    return()
endif()

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
