# We use GCC/Clang's function multi-versioning for AVX512
# support. This code will automatically dispatch to the
# AVX512 VPOPCNT algorithm if the CPU supports it and use
# the default (portable) algorithm otherwise.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

check_cxx_source_compiles("
    // GCC/Clang function multiversioning for AVX512 is not needed if
    // the user compiles with -mavx512f -mavx512vpopcntdq -mbmi2.
    // GCC/Clang function multiversioning generally causes a minor
    // overhead, hence we disable it if it is not needed.
    #if defined(__AVX512F__) && \
        defined(__AVX512VPOPCNTDQ__) && \
        defined(__BMI2__)
      Error: AVX512 BMI2 multiarch not needed!
    #endif

    #include <cpu_supports_avx512_bmi2.hpp>
    #include <immintrin.h>
    #include <stdint.h>

    class Sieve {
    public:
        uint64_t count_default(uint64_t* array, uint64_t stop_idx);
        __attribute__ ((target (\"avx512f,avx512vpopcntdq,bmi2\")))
        uint64_t count_avx512_bmi2(uint64_t* array, uint64_t stop_idx);
    };

    uint64_t Sieve::count_default(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t res = 0;
        for (uint64_t i = 0; i < stop_idx; i++)
            res += array[i];
        return res;
    }

    __attribute__ ((target (\"avx512f,avx512vpopcntdq,bmi2\")))
    uint64_t Sieve::count_avx512_bmi2(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t i = 0;
        __m512i vcnt = _mm512_setzero_si512();
        do
        {
          __mmask8 mask = (i + 8 < stop_idx) ? 0xff : (__mmask8) _bzhi_u64(0xff, stop_idx - i);
          __m512i vec = _mm512_maskz_loadu_epi64(mask , &array[i]);
          vec = _mm512_popcnt_epi64(vec);
          vcnt = _mm512_add_epi64(vcnt, vec);
          i += 8;
        }
        while (i < stop_idx);
        return _mm512_reduce_add_epi64(vcnt);
    }

    int main()
    {
        uint64_t array[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        uint64_t cnt = 0;
        Sieve sieve;

        if (cpu_supports_avx512_bmi2)
            cnt = sieve.count_avx512_bmi2(&array[0], 10);
        else
            cnt = sieve.count_default(&array[0], 10);

        return (cnt > 0) ? 0 : 1;
    }
" multiarch_avx512_vpopcnt)

if(multiarch_avx512_vpopcnt)
    set(ENABLE_MULTIARCH "ENABLE_MULTIARCH_AVX512_BMI2")
endif()

cmake_pop_check_state()
