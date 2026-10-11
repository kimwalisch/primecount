# Runtime dispatch selects AVX512 VPOPCNT when supported by the CPU
# and OS, and uses the default algorithm otherwise. The kernels use
# target attributes with GCC/Clang (including clang-cl), and noinline
# with native MSVC.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")

check_cxx_source_compiles("
    #include <cpu_arch_macros.hpp>

    // Runtime dispatch for AVX512 is not needed if
    // the user compiles with -mavx512f -mavx512bw -mavx512vl
    // -mavx512vpopcntdq -mpopcnt.
    // Runtime dispatch generally causes a minor
    // overhead, hence we disable it if it is not needed.
    #if defined(ENABLE_AVX512_VPOPCNT)
      Error: AVX512 multiarch not needed!
    #endif

    #include <macros.hpp>
    #include <src/arch/x86/cpuid.cpp>
    #include <immintrin.h>
    #include <stdint.h>

    class Sieve {
    public:
        uint64_t count_default(uint64_t* array, uint64_t stop_idx);
        PRIMECOUNT_TARGET_ENTRY(\"avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt\")
        uint64_t count_avx512(uint64_t* array, uint64_t stop_idx);
    };

    uint64_t Sieve::count_default(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t res = 0;
        for (uint64_t i = 0; i < stop_idx; i++)
            res += array[i];
        return res;
    }

    PRIMECOUNT_TARGET_ENTRY(\"avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt\")
    uint64_t Sieve::count_avx512(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t i = 0;
        __m512i vcnt = _mm512_setzero_si512();

        for (; i + 8 < stop_idx; i += 8)
        {
            __m512i vec = _mm512_loadu_epi64(&array[i]);
            vec = _mm512_popcnt_epi64(vec);
            vcnt = _mm512_add_epi64(vcnt, vec);
        }

        __mmask8 mask = (__mmask8) (0xff >> (i + 8 - stop_idx));
        __m512i vec = _mm512_maskz_loadu_epi64(mask , &array[i]);
        vec = _mm512_popcnt_epi64(vec);
        vcnt = _mm512_add_epi64(vcnt, vec);
        return _mm512_reduce_add_epi64(vcnt) + _mm_popcnt_u64(array[0]);
    }

    int main()
    {
        uint64_t array[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        uint64_t cnt = 0;
        Sieve sieve;

        if (primecount::has_avx512_vpopcnt())
            cnt = sieve.count_avx512(&array[0], 10);
        else
            cnt = sieve.count_default(&array[0], 10);

        return (cnt > 0) ? 0 : 1;
    }
" multiarch_avx512_vpopcnt)

if(multiarch_avx512_vpopcnt)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "ENABLE_MULTIARCH_AVX512_VPOPCNT")
endif()

cmake_pop_check_state()
