#ifndef __LIBKRCRAND_MATH_HPP__
#define __LIBKRCRAND_MATH_HPP__

#include <cstdint>
#include "platform.hpp"

namespace krcrand{
inline double type_depend_mull(double a, double b){
    return a*b;
}
double uniform01_exclude0(uint64_t a);
double unsafe_log(double x);
double unsafe_lambertw0(double x);
}
#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
namespace krcrand{
inline __m128d type_depend_mull(__m128d a, double b)
{
    return _mm_mul_pd(a, _mm_set1_pd(b));
}
__m128d uniform01_exclude0(__m128i a);
__m128d  unsafe_log(__m128d x);
}
#endif


#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
namespace krcrand{
inline __m256d type_depend_mull(__m256d a, double b)
{
    return _mm256_mul_pd(a, _mm256_set1_pd(b));
}
__m256d uniform01_exclude0(__m256i a);
__m256d  unsafe_log(__m256d x);
}
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
namespace krcrand{
inline __m512d type_depend_mull(__m512d a, double b)
{
    return _mm512_mul_pd(a, _mm512_set1_pd(b));
}
__m512d uniform01_exclude0(__m512i a);
__m512d  unsafe_log(__m512d x);
}
#endif
#endif