# Runtime dispatch selects AVX2 when supported by the CPU and OS,
# and uses the default algorithm otherwise. The kernels use target
# attributes with GCC/Clang (including clang-cl), and noinline with
# native MSVC.

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()
set(CMAKE_REQUIRED_INCLUDES
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")

check_cxx_source_compiles("
    #include <cpu_arch_macros.hpp>

    // Runtime dispatch for AVX2 is not needed if the compiler flags
    // already enable AVX2 and POPCNT, or the faster AVX512 VPOPCNT
    // algorithm. Runtime dispatch generally causes a minor
    // overhead, hence we disable it if it is not needed.
    #if defined(ENABLE_AVX2) || \
        defined(ENABLE_AVX512_VPOPCNT)
      Error: AVX2 multiarch not needed!
    #endif

    #include <macros.hpp>
    #include <src/arch/x86/cpuid.cpp>
    #include <immintrin.h>
    #include <stdint.h>

    uint32_t filter_default(const uint16_t* factor_table, uint32_t* m_indexes)
    {
        uint32_t m_count = 0;
        for (uint32_t m = 0; m < 8; m++)
        {
            m_indexes[m_count] = m;
            m_count += (factor_table[m] > 3);
        }
        return m_count;
    }

    PRIMECOUNT_TARGET_ENTRY(\"avx2,popcnt\")
    uint32_t filter_avx2(const uint16_t* factor_table, uint32_t* m_indexes)
    {
        const uint8_t offsets[8] = { 0, 2, 4, 6, 0, 0, 0, 0 };
        __m256i vec = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i*) factor_table));
        __m256i cmp = _mm256_cmpgt_epi32(vec, _mm256_set1_epi32(3));
        uint32_t mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp));
        __m256i offsets32 = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i*) offsets));
        __m256i indexes = _mm256_sub_epi32(_mm256_set1_epi32(7), offsets32);
        _mm256_storeu_si256((__m256i*) m_indexes, indexes);
        return _mm_popcnt_u32(mask);
    }

    int main()
    {
        uint16_t factor_table[8] = { 1, 5, 2, 7, 0, 9, 3, 4 };
        uint32_t m_indexes[8];
        uint32_t m_count;

        if (primecount::has_avx2())
            m_count = filter_avx2(factor_table, m_indexes);
        else
            m_count = filter_default(factor_table, m_indexes);

        return (m_count == 4) ? 0 : 1;
    }
" PRIMECOUNT_MULTIARCH_AVX2)

if(PRIMECOUNT_MULTIARCH_AVX2)
    list(APPEND PRIMECOUNT_COMPILE_DEFINITIONS "ENABLE_MULTIARCH_AVX2")
endif()

cmake_pop_check_state()
