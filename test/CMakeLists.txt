file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE ${PRIMECOUNT_COMPILE_DEFINITIONS})
    target_link_libraries(${binary_name} primecount::primecount primesieve::primesieve ${PRIMECOUNT_LINK_LIBRARIES})

    target_include_directories(${binary_name}
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/deleglise-rivat
        ${CMAKE_SOURCE_DIR}/src/gourdon
        ${CMAKE_SOURCE_DIR}/src/lmo)

    add_test(NAME ${binary_name} COMMAND ${binary_name})

    # Copy primecount.dll and primesieve.dll to test directory.
    # On Windows the DLLs must be in the same directory as the
    # binaries that depend on them.
    if (WIN32 AND NOT STATICALLY_LINK_LIBPRIMECOUNT)
        add_custom_command(TARGET ${binary_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:libprimecount>
                $<TARGET_FILE_DIR:${binary_name}>)

        add_custom_command(TARGET ${binary_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:libprimesieve>
                $<TARGET_FILE_DIR:${binary_name}>)
    endif()
endforeach()

# The order here is important, we want simpler implementations
# to be tested before more advanced/complex implementations.
# LMO is less complex than deleglise-rivat, deleglise-rivat is
# less complex than gourdon, and the C/C++ API depends on all
# these implementations.
add_subdirectory(lmo)
add_subdirectory(deleglise-rivat)
add_subdirectory(gourdon)
add_subdirectory(api)
