file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE "${HAVE_FLOAT128}" "${DISABLE_INT128}" "${ENABLE_DIV32}" "${ENABLE_ASSERT}" "${ENABLE_MULTIARCH}" "${ENABLE_INT128_OPENMP_PATCH}")
    target_link_libraries(${binary_name} primecount::primecount primesieve::primesieve "${LIB_OPENMP}" "${LIB_ATOMIC}")
    add_test(NAME ${binary_name} COMMAND ${binary_name})
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
