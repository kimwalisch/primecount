if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    return()
endif()

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" PRIMECOUNT_SYSTEM_PROCESSOR_LOWER)
string(TOLOWER "${CMAKE_CXX_COMPILER_ARCHITECTURE_ID}" PRIMECOUNT_COMPILER_ARCH_LOWER)

if(NOT PRIMECOUNT_SYSTEM_PROCESSOR_LOWER MATCHES "^(x86_64|amd64)$" AND
   NOT PRIMECOUNT_COMPILER_ARCH_LOWER MATCHES "^(x64|x86_64|amd64)$")
    return()
endif()

set(PRIMECOUNT_CODEGEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(PRIMECOUNT_CODEGEN_SRC "${CMAKE_CURRENT_SOURCE_DIR}/conditional_move.cpp")
set(PRIMECOUNT_CODEGEN_HEADERS
    "${PROJECT_SOURCE_DIR}/include/macros.hpp"
    "${PROJECT_SOURCE_DIR}/include/popcnt.hpp")
set(PRIMECOUNT_CODEGEN_OBJ "${PRIMECOUNT_CODEGEN_DIR}/conditional_move.obj")
set(PRIMECOUNT_CODEGEN_ASM "${PRIMECOUNT_CODEGEN_DIR}/conditional_move.s")
set(PRIMECOUNT_CODEGEN_FLAGS "${CMAKE_CXX_FLAGS}")
set(PRIMECOUNT_CODEGEN_CHECKER
    "${CMAKE_CURRENT_SOURCE_DIR}/check_conditional_move_codegen.cmake")
separate_arguments(PRIMECOUNT_CODEGEN_FLAGS NATIVE_COMMAND "${PRIMECOUNT_CODEGEN_FLAGS}")

if(MSVC)
    set(PRIMECOUNT_CODEGEN_ASM "${PRIMECOUNT_CODEGEN_DIR}/conditional_move.asm")
    add_custom_command(
        OUTPUT "${PRIMECOUNT_CODEGEN_ASM}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${PRIMECOUNT_CODEGEN_DIR}"
        COMMAND "${CMAKE_CXX_COMPILER}"
                /nologo
                /std:c++17
                /O2
                /FAs
                /c
                ${PRIMECOUNT_CODEGEN_FLAGS}
                "/I${PROJECT_SOURCE_DIR}/include"
                "/Fo${PRIMECOUNT_CODEGEN_OBJ}"
                "/Fa${PRIMECOUNT_CODEGEN_ASM}"
                "${PRIMECOUNT_CODEGEN_SRC}"
        DEPENDS "${PRIMECOUNT_CODEGEN_SRC}" ${PRIMECOUNT_CODEGEN_HEADERS}
        COMMENT "Generating assembly for check_conditional_move_codegen"
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
        COMMENT "Generating assembly for check_conditional_move_codegen"
        VERBATIM)
endif()

add_custom_target(check_conditional_move_codegen
    COMMAND ${CMAKE_COMMAND}
            "-DASM_FILE=${PRIMECOUNT_CODEGEN_ASM}"
            -P "${PRIMECOUNT_CODEGEN_CHECKER}"
    DEPENDS "${PRIMECOUNT_CODEGEN_ASM}"
            "${PRIMECOUNT_CODEGEN_CHECKER}"
    VERBATIM)

add_test(
    NAME check_conditional_move_codegen
    COMMAND ${CMAKE_COMMAND}
            --build "${PROJECT_BINARY_DIR}"
            --target check_conditional_move_codegen
            --config "$<CONFIG>")
