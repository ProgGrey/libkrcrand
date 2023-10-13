#include "tools.hpp"
#include <cstring>

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#endif

#ifdef LIBKRCRAND_ENABLE_SSE41
#include <tmmintrin.h>
#include <smmintrin.h>
#endif
void split_mix64_fill(uint64_t seed, uint64_t *buf, unsigned int size)
{
    for(unsigned int k = 0; k < size; k++){
        uint64_t val = (seed += 0x9e3779b97f4a7c15);
        val = (val ^ (val >> 30)) * 0xbf58476d1ce4e5b9;
        val = (val ^ (val >> 27)) * 0x94d049bb133111eb;
        buf[k] =  val ^ (val >> 31);
    }
}

uint64_t rotl(const uint64_t x, uint_fast8_t k) 
{
	return (x << k) | (x >> (64 - k));
}

#ifdef LIBKRCRAND_ENABLE_SSE2
__m128i rotl_sse2(const __m128i x, uint_fast8_t k) 
{
    auto a = _mm_slli_epi64(x, k);
    auto b = _mm_srli_epi64(x, 64 - k);
	return _mm_or_si128(a, b);
}

__m128i mull_u64_sse2(const __m128i &a, const __m128i &b)
{
    /*
    __m128i a1 = a;
    __m128i b1 = b;
    __m128i t1 = _mm_mul_epu32(a1, b1);
    __m128i m23 = _mm_mullo_epi32(a, _mm_shuffle_epi32(b, 0b10110001));
    __m128i t2 = _mm_hadd_epi32(m23, _mm_set1_epi64x(0));
    t2 = _mm_shuffle_epi32(t2, 0b01110010);
    return _mm_add_epi64(t1, t2);
    //*/
    
    // It's faster that direct solution via 3 multiplication
    uint64_t ab[2];
    uint64_t bb[2];
    memcpy(ab, &a, sizeof(uint64_t)*2);
    memcpy(bb, &b, sizeof(uint64_t)*2);
    ab[0] *= bb[0];
    ab[1] *= bb[1];
    __m128i res;
    memcpy(&res, ab, sizeof(uint64_t)*2);
    return res;


    /*
    //__m128i a1 = _mm_and_si128(a, _mm_set1_epi64x(0xFFFFFFFF));
    //__m128i b1 = _mm_and_si128(b, _mm_set1_epi64x(0xFFFFFFFF));
    __m128i a1 = a;
    __m128i b1 = b;
    __m128i a2 = _mm_srli_epi64(a, 32);
    __m128i b2 = _mm_srli_epi64(b, 32);
    __m128i t1 = _mm_mul_epu32(a1, b1);
    __m128i t2 = _mm_slli_epi64(_mm_mul_epu32(a1, b2), 32);
    __m128i t3 = _mm_slli_epi64(_mm_mul_epu32(a2, b1), 32);
    return _mm_add_epi64(_mm_add_epi64(t1, t2), t3);
    //*/
}
#endif

#ifdef LIBKRCRAND_ENABLE_AVX2
__m256i rotl_avx2(const __m256i x, uint_fast8_t k) 
{
    auto a = _mm256_slli_epi64(x, k);
    auto b = _mm256_srli_epi64(x, 64 - k);
	return _mm256_or_si256(a, b);
}
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
__m512i rotl_avx512f(const __m512i x, uint_fast8_t k) 
{
    auto a = _mm512_slli_epi64(x, k);
    auto b = _mm512_srli_epi64(x, 64 - k);
	return _mm512_or_si512(a, b);
}
#endif