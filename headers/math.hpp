#ifndef __LIBKRCRAND_MATH_HPP__
#define __LIBKRCRAND_MATH_HPP__

#include <cstdint>
#include "platform.hpp"

namespace krcrand{

const double M_E_m1 = 0.367879441171442321595523770162;//e^-1
const double M_E_m2 = 0.135335283236612691893999494973;//e^-2

// solve the equation a*x^3+b*x^2+c*x+d=0
// @return number of real roots.
int p3_solve(double a, double b, double c, double d, double &x1, double &x2, double &x3);
int p3_solve_cardano(long double const &a, long double const &b, long double const &c, long double const &d, double &x1, double &x2, double &x3);
int p3_solve_modern(long double const &a, long double const &b, long double const &c, long double const &d, double &x1, double &x2, double &x3);

// solve the equation a*x^2+b*x+c=0
// @return number of real roots.
int p2_solve(double a, double b, double c, double &x1, double &x2);

// solve the equation 0.5a*x^2+b*x+c=0
double unsafe_p2_half_solve(double a, double b, double c, double xl, double xr);

double signum(double x);

inline double type_depend_mull(double a, double b){
    return a*b;
}

double uniform01_exclude01(uint64_t a);
double uniform01_exclude0(uint64_t a);
double uniform01(uint64_t a);
double unsafe_log(double x);
double unsafe_lambertw0(double x);
double unsafe_lambertw1(double x);
}
#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
namespace krcrand{
inline __m128d type_depend_mull(__m128d a, double b)
{
    return _mm_mul_pd(a, _mm_set1_pd(b));
}
__m128d uniform01(__m128i a);
__m128d uniform01_exclude0(__m128i a);
__m128d uniform01_exclude01(__m128i a);
__m128d  unsafe_log(__m128d x);
__m128d unsafe_p2_half_solve(__m128d a, __m128d b, __m128d c, __m128d xl, __m128d xr);
}
#endif


#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
namespace krcrand{
inline __m256d type_depend_mull(__m256d a, double b)
{
    return _mm256_mul_pd(a, _mm256_set1_pd(b));
}
__m256d uniform01(__m256i a);
__m256d uniform01_exclude0(__m256i a);
__m256d uniform01_exclude01(__m256i a);
__m256d  unsafe_log(__m256d x);
__m256d unsafe_p2_half_solve(__m256d a, __m256d b, __m256d c, __m256d xl, __m256d xr);
}
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
namespace krcrand{
inline __m512d type_depend_mull(__m512d a, double b)
{
    return _mm512_mul_pd(a, _mm512_set1_pd(b));
}
__m512d uniform01(__m512i a);
__m512d uniform01_exclude0(__m512i a);
__m512d uniform01_exclude01(__m512i a);
__m512d unsafe_log(__m512d x);
__m512d unsafe_p2_half_solve(__m512d a, __m512d b, __m512d c, __m512d xl, __m512d xr);
}
#endif
#endif