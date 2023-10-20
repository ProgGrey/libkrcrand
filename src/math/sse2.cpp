#include "../../headers/platform.hpp"
#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#include "../math.hpp"

#define d_type __m128d
#define i_type __m128i

#define mull_d(a, b) _mm_mul_pd(a, b)
#define div_d(a,b) _mm_div_pd(a, b)
#define add_d(a, b) _mm_add_pd(a, b)
#define sub_d(a, b) _mm_sub_pd(a, b)
#ifdef LIBKRCRAND_ENABLE_FMA3
#include <immintrin.h>
#define fma_d(a,b,c) _mm_fmadd_pd(a,b,c)
#else
#define fma_d(a,b,c) _mm_add_pd(_mm_mul_pd(a,b), c)
#endif

#define to_i(a) _mm_castpd_si128(a)
#define to_d(a) _mm_castsi128_pd(a)

#define shr_i64(a,b) _mm_srli_epi64(a, b)
#define shl_i64(a,b) _mm_slli_epi64(a, b)

#define set1_d(a) _mm_set1_pd(a)
#define set1_i(a) _mm_set1_epi64x(a)

#define le_d(a,b) _mm_cmple_pd(a, b)
#define lt_d(a,b) _mm_cmplt_pd(a, b)
#define notand_d(a,b) _mm_andnot_pd(a, b)
#define and_d(a,b) _mm_and_pd(a, b)
#define or_d(a,b) _mm_or_pd(a,b)
#define and_i(a,b) _mm_and_si128(a,b)
#define or_i(a,b) _mm_or_si128(a,b)

#include "amd64.hpp"

horner_function(horner_sse2)
horner1_function(horner1_sse2)
log_function(unsafe_log_sse2, horner_sse2, horner1_sse2)

#endif