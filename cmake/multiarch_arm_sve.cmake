# We use GCC/Clang's function multi-versioning for ARM SVE
# support. This code will automatically dispatch to the
# ARM SVE algorithm if the CPU supports it and use the default
# (portable) algorithm otherwise.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}")

check_cxx_source_compiles("
    // GCC/Clang function multiversioning for ARM SVE is not needed
    // if the user compiles with -march=armv8-a+sve. GCC/Clang
    // function multiversioning generally causes a minor overhead,
    // hence we disable it if it is not needed.
    #if defined(__ARM_FEATURE_SVE) && \
        __has_include(<arm_sve.h>)
      Error: ARM SVE multiarch not needed!
    #endif

    #include <src/arch/arm/sve.cpp>
    #include <arm_sve.h>
    #include <stdint.h>

    class Sieve {
    public:
        uint64_t count_default(uint64_t* array, uint64_t stop_idx);
        __attribute__ ((target (\"arch=armv8-a+sve\")))
        uint64_t count_arm_sve(uint64_t* array, uint64_t stop_idx);
    };

    uint64_t Sieve::count_default(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t res = 0;
        for (uint64_t i = 0; i < stop_idx; i++)
            res += array[i];
        return res;
    }

    __attribute__ ((target (\"arch=armv8-a+sve\")))
    uint64_t Sieve::count_arm_sve(uint64_t* array, uint64_t stop_idx)
    {
        svuint64_t vec = svinsr_n_u64(svdup_u64(array[0]), array[1]);
        svuint64_t vcnt = svcnt_u64_z(svwhilelt_b64(0, 2), vec);
        uint64_t i = 2;

        for (; i + svcntd() < stop_idx; i += svcntd())
        {
            vec = svld1_u64(svptrue_b64(), &array[i]);
            vec = svcnt_u64_x(svptrue_b64(), vec);
            vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
        }

        svbool_t pg = svwhilelt_b64(i, stop_idx);
        vec = svld1_u64(pg, &array[i]);
        vec = svcnt_u64_z(pg, vec);
        vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
        return svaddv_u64(svptrue_b64(), vcnt);
    }

    int main()
    {
        uint64_t array[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        uint64_t cnt = 0;
        Sieve sieve;

        if (primecount::has_arm_sve())
            cnt = sieve.count_arm_sve(&array[0], 10);
        else
            cnt = sieve.count_default(&array[0], 10);

        return (cnt > 0) ? 0 : 1;
    }
" multiarch_arm_sve)

if(multiarch_arm_sve)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "ENABLE_MULTIARCH_ARM_SVE")
endif()

cmake_pop_check_state()
