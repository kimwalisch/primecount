include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
    #include <limits>
    #include <type_traits>
    int main() {
        static_assert(std::numeric_limits<unsigned int>::min() == 0, \"\");
        static_assert(std::is_integral<int>::value, \"\");
        auto lambda = [](auto x) { return x + 1; };
        return (int) lambda(0) - 1;
    }" compiler_supports_cpp14)
