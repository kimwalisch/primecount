# Apple Silicon CPUs first released in 2020 have very fast integer
# division instructions. For Apple Silicon CPUs on Apple OSes
# we disable libdivide to get the best performance.

include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
    #if defined(__APPLE__) && defined(__aarch64__)
        int main() { return 0; }
    #else
        Error: not Apple ARM64
    #endif
    " Apple_ARM64)
