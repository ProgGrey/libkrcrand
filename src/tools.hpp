#ifndef __LIBKRCRAND_TOOLS_HPP__
#define __LIBKRCRAND_TOOLS_HPP__

#include <cstdint>
#include "../headers/platform.hpp"

void split_mix64_fill(uint64_t seed, uint64_t *buf, unsigned int  size);
uint64_t rotl(const uint64_t x, uint_fast8_t k);

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
__m128i rotl_sse2(const __m128i x, uint_fast8_t k) ;
__m128i mull_u64_sse2(const __m128i &a, const __m128i &b);
#endif

#endif

