#include "tools.hpp"


#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
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
__m128i rotl_SSE2(const __m128i x, uint_fast8_t k) 
{
    auto a = _mm_slli_epi64(x, k);
    auto b = _mm_srli_epi64(x, 64 - k);
	return _mm_or_si128(a, b);
}
#endif