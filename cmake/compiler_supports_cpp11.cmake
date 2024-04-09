include(CheckCXXSourceCompiles)

# Check if compiler supports C++11 or later
check_cxx_source_compiles("
    #include <limits>
    #include <type_traits>
    int main() {
        static_assert(std::numeric_limits<unsigned int>::min() == 0, \"\");
        static_assert(std::is_integral<int>::value, \"\");
        return 0;
    }" compiler_supports_cpp11)
