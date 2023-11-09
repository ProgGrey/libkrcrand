#include "../../headers/platform.hpp"
#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
#include "../../headers/math.hpp"

using namespace krcrand;

#define d_type __m256d
#define i_type __m256i

#define mull_d(a, b) _mm256_mul_pd(a, b)
#define div_d(a,b) _mm256_div_pd(a, b)
#define add_d(a, b) _mm256_add_pd(a, b)
#define sub_d(a, b) _mm256_sub_pd(a, b)
#ifdef LIBKRCRAND_ENABLE_FMA3
#define fma_d(a,b,c) _mm256_fmadd_pd(a,b,c)
#else
#define fma_d(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b), c)
#endif

#define to_i(a) _mm256_castpd_si256(a)
#define to_d(a) _mm256_castsi256_pd(a)

#define shr_i64(a,b) _mm256_srli_epi64(a, b)
#define shl_i64(a,b) _mm256_slli_epi64(a, b)

#define set1_d(a) _mm256_set1_pd(a)
#define set1_i(a) _mm256_set1_epi64x(a)

#define le_d(a,b) _mm256_cmp_pd (a, b, _CMP_LE_OS)
#define lt_d(a,b) _mm256_cmp_pd (a, b, _CMP_LT_OS)
#define notand_d(a,b) _mm256_andnot_pd(a, b)
#define and_d(a,b) _mm256_and_pd(a, b)
#define or_d(a,b) _mm256_or_pd(a,b)
#define and_i(a,b) _mm256_and_si256(a,b)
#define or_i(a,b) _mm256_or_si256(a,b)


#include "amd64.hpp"

#if defined(LIBKRCRAND_ENABLE_AVX512VL) && defined(LIBKRCRAND_ENABLE_AVX512DQ)
#define  conv_u_to_d(x) _mm256_cvtepu64_pd(x)
#else
//Based on solution from https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
__m256d u64_to_d_avx2(const __m256i x)
{
    __m256i magic_i_lo = _mm256_set1_epi64x(0x4330000000000000);// 2^52  encoded as floating-point 
    __m256i magic_i_hi32 = _mm256_set1_epi64x(0x4530000000000000); // 2^84 encoded as floating-point 
    __m256i magic_i_all = _mm256_set1_epi64x(0x4530000000100000); // 2^84 + 2^52 encoded as floating-point
    __m256d magic_d_all = _mm256_castsi256_pd(magic_i_all);

    __m256i x_lo = _mm256_blend_epi32(magic_i_lo, x, 0b01010101);// Blend the 32 lowest significant bits of x with magic_int_lo
    __m256i x_hi = _mm256_srli_epi64(x, 32); // Extract the 32 most significant bits of x
            x_hi = _mm256_xor_si256(x_hi, magic_i_hi32);// Blend v_hi with 0x45300000
    __m256d x_hi_dbl = _mm256_sub_pd(_mm256_castsi256_pd(x_hi), magic_d_all); // Compute in double precision:
    __m256d result = _mm256_add_pd(x_hi_dbl, _mm256_castsi256_pd(x_lo));// (x_hi - magic_d_all) + x_lo  Do not assume associativity of floating point addition !!
    return result;
}
#define  conv_u_to_d(x) u64_to_d_avx2(x)
#endif

horner_function(horner_avx2)
horner1_function(horner1_avx2)

namespace krcrand{
uniform01_exclude0_function(uniform01_exclude0)
log_function(unsafe_log, horner_avx2, horner1_avx2)
}

#endif