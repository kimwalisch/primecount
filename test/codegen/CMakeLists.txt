file(GLOB codegen_tests CONFIGURE_DEPENDS "*.cmake")

set(PRIMECOUNT_CONFIGURE_CODEGEN_TEST ON)

foreach(codegen_test IN LISTS codegen_tests)
    include("${codegen_test}")
endforeach()

unset(PRIMECOUNT_CONFIGURE_CODEGEN_TEST)
