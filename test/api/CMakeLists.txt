enable_language(C)
file(GLOB files "*.cpp" "*.c")

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
endforeach()
