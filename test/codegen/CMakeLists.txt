file(GLOB codegen_sources CONFIGURE_DEPENDS "*.cpp")

set(PRIMECOUNT_CONFIGURE_CODEGEN_TEST ON)

foreach(codegen_source IN LISTS codegen_sources)
    get_filename_component(codegen_name "${codegen_source}" NAME_WE)
    set(codegen_test "${CMAKE_CURRENT_SOURCE_DIR}/${codegen_name}.cmake")

    if(NOT EXISTS "${codegen_test}")
        message(FATAL_ERROR "Missing codegen module: ${codegen_test}")
    endif()

    include("${codegen_test}")
endforeach()

unset(PRIMECOUNT_CONFIGURE_CODEGEN_TEST)
