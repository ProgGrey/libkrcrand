#include "../../headers/platform.hpp"
#ifdef LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
#include "../../headers/math.hpp"

using namespace krcrand;

#define d_type __m512d
#define i_type __m512i

#define mull_d(a, b) _mm512_mul_pd(a, b)
#define div_d(a,b) _mm512_div_pd(a, b)
#define add_d(a, b) _mm512_add_pd(a, b)
#define sub_d(a, b) _mm512_sub_pd(a, b)
#define fma_d(a,b,c) _mm512_fmadd_pd(a,b,c)

#define to_i(a) _mm512_castpd_si512(a)
#define to_d(a) _mm512_castsi512_pd(a)

#define shr_i64(a,b) _mm512_srli_epi64(a, b)
#define shl_i64(a,b) _mm512_slli_epi64(a, b)

#define set1_d(a) _mm512_set1_pd(a)
#define set1_i(a) _mm512_set1_epi64(a)

#define le_d(a,b) _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define lt_d(a,b) _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define notand_d(a,b) _mm512_andnot_pd(a, b)
#define and_d(a,b) _mm512_and_pd(a, b)
#define or_d(a,b) _mm512_or_pd(a,b)
#define and_i(a,b) _mm512_and_si512(a,b)
#define or_i(a,b) _mm512_or_si512(a,b)

#define mov_mask_d(flag, if_true, if_false) _mm512_mask_mov_pd (if_false, flag, if_true)

#define __AVX512_BRANCH__
#include "amd64.hpp"

#ifdef LIBKRCRAND_ENABLE_AVX512DQ
#define  conv_u_to_d(x) _mm512_cvtepu64_pd(x)
#else
__m512d u64_to_d_avx512(__m512i x)
{
    alignas(64) uint64_t val[8];
    _mm512_store_epi64(val, x);
    alignas(64) double res[8];
    for(unsigned int k =0; k < 8; k++){
        res[k] = val[k];
    }
    return _mm512_load_pd(res);
    
}
#define  conv_u_to_d(x) u64_to_d_avx512(x)
#endif

horner_function(horner_avx512)
horner1_function(horner1_avx512)
namespace krcrand{
uniform01_exclude0_function(uniform01_exclude0)
log_function(unsafe_log, horner_avx512, horner1_avx512)
}
#endif