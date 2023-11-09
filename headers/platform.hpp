#ifndef __LIBKRCRAND_PLATFORM_HPP__
#define __LIBKRCRAND_PLATFORM_HPP__

namespace krcrand{

#ifdef __SSE2__
#define LIBKRCRAND_ENABLE_SSE2 true
#endif
#ifdef __SSE4_1__
#define LIBKRCRAND_ENABLE_SSE4_1 true
#endif
#ifdef __AVX2__
#define LIBKRCRAND_ENABLE_AVX2 true
#endif
#ifdef __AVX512F__ 
#define LIBKRCRAND_ENABLE_AVX512F true
#endif
#ifdef __AVX512VL__
#define LIBKRCRAND_ENABLE_AVX512VL true
#endif
#ifdef __AVX512DQ__
#define LIBKRCRAND_ENABLE_AVX512DQ true
#endif

#ifdef __FMA__
#define LIBKRCRAND_ENABLE_FMA3 true
#endif


#ifdef LIBKRCRAND_ENABLE_AVX512F
#define DECL_KRCRAND_ALIGN alignas(64)
// Buffer size for entropy generators. Must be lesser than uint16_t maximum value (65536) and can be divided by 16. Different values below for better performance
#define Generator_Buff_Size 128

#elif defined(LIBKRCRAND_ENABLE_AVX2)
#define DECL_KRCRAND_ALIGN alignas(32)
#define Generator_Buff_Size 32

#elif defined(LIBKRCRAND_ENABLE_SSE2)
#define Generator_Buff_Size 16
#define DECL_KRCRAND_ALIGN alignas(16)

#else
#define DECL_KRCRAND_ALIGN 
#endif

#define COMPAB_COUNT_UNIV 8
#define COMPAB_COUNT_SSE2 4
#define COMPAB_COUNT_AVX2 2
#define COMPAB_COUNT_AVX512F 1

}

#endif