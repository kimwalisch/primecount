file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE ${PRIMECOUNT_COMPILE_DEFINITIONS})
    target_link_libraries(${binary_name} primecount::primecount primesieve::primesieve "${LIB_OPENMP}" "${LIB_ATOMIC}")
    add_test(NAME ${binary_name} COMMAND ${binary_name})
endforeach()
