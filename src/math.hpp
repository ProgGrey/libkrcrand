#ifndef __LIBKRCRAND_TOOLS_HPP__
#define __LIBKRCRAND_TOOLS_HPP__

#include "../headers/platform.hpp"

double unsafe_log(double x);

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
__m128d  unsafe_log_sse2(__m128d x);
#endif


#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
__m256d  unsafe_log_avx2(__m256d x);
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
__m512d  unsafe_log_avx512(__m512d x);
#endif


#endif