file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE "${DISABLE_POPCNT}")
    target_link_libraries(${binary_name} libprimecount primesieve::primesieve "${LIB_OPENMP}")
    add_test(NAME ${binary_name} COMMAND ${binary_name})
endforeach()
