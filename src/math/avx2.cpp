#include "../../headers/platform.hpp"
#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
#include "../math.hpp"

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

horner_function(horner_avx2)
horner1_function(horner1_avx2)
log_function(unsafe_log_avx2, horner_avx2, horner1_avx2)
#endif