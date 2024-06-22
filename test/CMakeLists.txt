file(GLOB files "*.cpp")

foreach(file ${files})
    get_filename_component(binary_name ${file} NAME_WE)
    add_executable(${binary_name} ${file})
    target_compile_definitions(${binary_name} PRIVATE ${PRIMECOUNT_COMPILE_DEFINITIONS})
    target_link_libraries(${binary_name} primecount::primecount primesieve::primesieve ${PRIMECOUNT_LINK_LIBRARIES})
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
