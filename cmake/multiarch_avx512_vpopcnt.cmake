include(CheckCXXSourceCompiles)

# We use GCC/Clang's function multi-versioning for AVX512
# support. This code will automatically dispatch to the
# AVX512 VPOPCNT algorithm if the CPU supports it and use
# the default (portable) algorithm otherwise.
check_cxx_source_compiles("
    #include <immintrin.h>
    #include <stdint.h>
    class Sieve {
        public:
        __attribute__ ((target (\"default\")))
        uint64_t count(uint64_t* array, uint64_t stop_idx);
        __attribute__ ((target (\"avx512f,avx512vpopcntdq\")))
        uint64_t count(uint64_t* array, uint64_t stop_idx);
    };
    __attribute__ ((target (\"default\")))
    uint64_t Sieve::count(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t res = 0;
        for (uint64_t i = 0; i < stop_idx; i++)
            res += array[i];
        return res;
    }
    __attribute__ ((target (\"avx512f,avx512vpopcntdq\")))
    uint64_t Sieve::count(uint64_t* array, uint64_t stop_idx)
    {
        uint64_t i = 0;
        __m512i vcnt = _mm512_setzero_si512();
        for (; i + 8 < stop_idx; i += 8)
        {
            __m512i vec = _mm512_loadu_epi64(&array[i]);
            vec = _mm512_popcnt_epi64(vec);
            vcnt = _mm512_add_epi64(vcnt, vec);
        }
        __mmask8 mask = 0xff >> (i + 8 - stop_idx);
        __m512i vec = _mm512_maskz_loadu_epi64(mask, &array[i]);
        vec = _mm512_popcnt_epi64(vec);
        vcnt = _mm512_add_epi64(vcnt, vec);
        return _mm512_reduce_add_epi64(vcnt);
    }
    int main()
    {
        uint64_t array[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        Sieve sieve;
        sieve.count(&array[0], 10);
        return 0;
    }
" multiarch_avx512_vpopcnt)
