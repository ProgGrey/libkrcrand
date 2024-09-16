#include "../../headers/platform.hpp"
#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#include "../../headers/math.hpp"

using namespace krcrand;

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
#define and_cond(a, b) _mm_and_pd(a, b)

#define notand_d(a,b) _mm_andnot_pd(a, b)
#define and_d(a,b) _mm_and_pd(a, b)
#define or_d(a,b) _mm_or_pd(a,b)
#define and_i(a,b) _mm_and_si128(a,b)
#define or_i(a,b) _mm_or_si128(a,b)
#define xor_d(a,b) _mm_xor_pd(a,b)

#define sqrt_d(a) _mm_sqrt_pd(a)

/*
#include <iostream>

void prnt_dc_sse(__m128d val)
{
    __m128i vi = to_i(val);
    alignas(16) uint64_t buf[2];
    _mm_store_si128(reinterpret_cast<__m128i*>(buf), vi);
    std::cout << std::hex;
    std::cout << "0x" << buf[0] << " 0x" << buf[1] << '\n';
    std::cout << std::dec;
}

#define prnt_dc(a) prnt_dc_sse(a)
//*/

#include "amd64.hpp"

#if defined(LIBKRCRAND_ENABLE_AVX512VL) && defined(LIBKRCRAND_ENABLE_AVX512DQ)
#define  conv_u_to_d(x) _mm_cvtepu64_pd(x)
#elif defined(LIBKRCRAND_ENABLE_SSE4_1)
#include <smmintrin.h>
// Based on solution from https://stackoverflow.com/questions/28939652/how-to-detect-sse-sse2-avx-avx2-avx-512-avx-128-fma-kcvi-availability-at-compile
__m128d u64_to_d_sse2(__m128i x){
    __m128i xH = _mm_srli_epi64(x, 32);
    xH = _mm_or_si128(xH, _mm_castpd_si128(_mm_set1_pd(19342813113834066795298816.)));          //  2^84
    __m128i xL = _mm_blend_epi16(x, _mm_castpd_si128(_mm_set1_pd(0x0010000000000000)), 0xcc);   //  2^52
    __m128d f = _mm_sub_pd(_mm_castsi128_pd(xH), _mm_set1_pd(19342813118337666422669312.));     //  2^84 + 2^52
    return _mm_add_pd(f, _mm_castsi128_pd(xL));
}
#define  conv_u_to_d(x) u64_to_d_sse2(x)
#else
__m128d u64_to_d_sse2(__m128i x){
    alignas(16) uint64_t val[2];
    alignas(16) double res[2];
    _mm_store_si128(reinterpret_cast<__m128i*>(val), x);
    res[0] = val[0];
    res[1] = val[1];
    return _mm_load_pd(res);
}
#define  conv_u_to_d(x) u64_to_d_sse2(x)
#endif

#ifdef LIBKRCRAND_ENABLE_SSE4_1
#include <smmintrin.h>
#define testc(a,b) _mm_testc_si128(_mm_castpd_si128(a),_mm_castpd_si128(b))
#define test_all_ones(x) _mm_test_all_ones(_mm_castpd_si128(x))
#else
int testc_sse2_slow(__m128i a, __m128i b){
    __m128i res = _mm_andnot_si128 (a, b);
    alignas(16) uint64_t res_arr[2];
    _mm_store_si128(reinterpret_cast<__m128i*>(res_arr), res);
    if((res_arr[0] == 0) && (res_arr[1] == 0)){
        return 1;
    } else{
        return 0;
    }
}
#define testc(a,b) testc_sse2_slow(_mm_castpd_si128(a), _mm_castpd_si128(b))
#define test_all_ones(x) testc_sse2_slow(_mm_castpd_si128(x), _mm_set1_epi64x(0xFFFFFFFFFFFFFFFF))
#endif


horner_function(horner_sse2)
horner1_function(horner1_sse2)
namespace krcrand{
uniform01_function(uniform01)
uniform01_exclude0_function(uniform01_exclude0)
uniform01_exclude01_function(uniform01_exclude01)
log_function(unsafe_log, horner_sse2, horner1_sse2)

p2_half_function(unsafe_p2_half_solve)
}
#endif