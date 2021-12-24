file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE "${DISABLE_POPCNT}" "${DISABLE_INT128}" "${ENABLE_DIV32}")
    target_link_libraries(${binary_name} libprimecount primesieve::primesieve "${LIB_OPENMP}" "${LIB_ATOMIC}")
    add_test(NAME ${binary_name} COMMAND ${binary_name})
endforeach()
