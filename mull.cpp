#include <cstdint>
#include <iostream>
#include <emmintrin.h>

using namespace std;

uint64_t part_mull(uint64_t a, uint64_t b)
{
    uint64_t a1 = a & 0xFFFFFFFF;    
    uint64_t b1 = b & 0xFFFFFFFF;
    uint64_t a2 = (a >> 32);
    uint64_t b2 = (b >> 32);
    uint64_t t1 = a1*b1;
    uint64_t t2 = a1*b2 << 32;
    uint64_t t3 = a2*b1 << 32;
    return t1+t2+t3;
}

__m128i mull_u64_sse2(const __m128i &a, const __m128i &b)
{
    __m128i a1 = _mm_and_si128(a, _mm_set1_epi64x(0xFFFFFFFF));
    __m128i b1 = _mm_and_si128(b, _mm_set1_epi64x(0xFFFFFFFF));
    __m128i a2 = _mm_srli_epi64(a, 32);
    __m128i b2 = _mm_srli_epi64(b, 32);
    __m128i t1 = _mm_mul_epu32(a1, b1);
    __m128i t2 = _mm_slli_epi64(_mm_mul_epu32(a1, b2), 32);
    __m128i t3 = _mm_slli_epi64(_mm_mul_epu32(a2, b1), 32);
    return _mm_add_epi64(_mm_add_epi64(t1, t2), t3);
}

int main()
{
    uint64_t a = 0xFFFFFFFFFFFFFFFF;
    uint64_t b = 0xFFFFFFFFFFFFFFFF;
    __m128i 
    cout << a*b << '\t' << part_mull(a,b) << '\t' << a*b-part_mull(a,b)<< endl;
    return 0;
}
