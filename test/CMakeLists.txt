file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_options(${binary_name} PRIVATE "${OpenMP_FLAGS}")
    target_link_libraries(${binary_name} libprimecount libprimesieve "${OpenMP_FLAGS}")
    add_test(NAME ${binary_name} COMMAND ${binary_name})
endforeach()
