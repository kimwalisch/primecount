file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_link_libraries(${binary_name} libprimecount libprimesieve)
    add_test(NAME ${binary_name} COMMAND ${binary_name})
endforeach()

target_compile_options(factor_table PRIVATE "${OpenMP_FLAGS}")

if(NOT MSVC_VERSION)
    target_link_libraries(factor_table "${OpenMP_FLAGS}")
endif()
