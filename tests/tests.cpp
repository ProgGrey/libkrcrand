#include "../headers/libkrcrand.hpp"
#include "../src/math.hpp"

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cmath>

using namespace std;
using namespace krcrand;


BOOST_AUTO_TEST_CASE(Xoshiro256mmGenerators)
{
    Xoshiro256mmState st;
    st.seed(1);
    const uint64_t st1[] = {0x910a2dec89025cc1, 0xbeeb8da1658eec67, 0xf893a2eefb32555e, 0x71c18690ee42c90b};
    const uint64_t st2[] = {0x53d630076a137ded, 0xed07f666882edfc6, 0x963ec9617b0bdbd3, 0x84b96906e4b2569a};    
    for(unsigned int k = 0; k < st.State_Size; k++){
        BOOST_CHECK(st.raw()[k] == st1[k]);
    }
    auto tmp = st.jump();
    for(unsigned int k = 0; k < st.State_Size; k++){
        BOOST_CHECK(st.raw()[k] == st1[k]);
    }
    for(unsigned int k = 0; k < tmp.State_Size; k++){
        BOOST_CHECK(tmp.raw()[k] == st2[k]);
    }
    Xoshiro256mmUniversalStable gen_u(st);
#ifdef LIBKRCRAND_ENABLE_SSE2
    Xoshiro256mmSSE2stable gen_sse(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    Xoshiro256mmAVX2stable gen_avx(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    Xoshiro256mmAVX512Fstable gen_avx512(st);
#endif
    for(unsigned int k = 0; k < 10000; k++){
        uint64_t val = gen_u();
#ifdef LIBKRCRAND_ENABLE_SSE2
        BOOST_CHECK(val == gen_sse());
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
        BOOST_CHECK(val == gen_avx());
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
        BOOST_CHECK(val == gen_avx512());
#endif
    }
    cout << endl;
}


BOOST_AUTO_TEST_CASE(math_tests)
{
    double values[8] = {0.1, 2, 3.5, 1, 0.000001, 12, 7643921, 2.2250738585072009e-318};
    double std_res[8];
    double uni_res[8];
    for(unsigned int k = 0; k < 8; k++){
        std_res[k] = log(values[k]);
        uni_res[k] = unsafe_log(values[k]);
        BOOST_CHECK(std_res[k] - uni_res[k] < 1e-15);
    }
#ifdef LIBKRCRAND_ENABLE_SSE2
    double sse2_res[8];
    for(unsigned int k = 0; k < 8; k+=2){
        __m128d val = _mm_loadu_pd(values + k);
        val = unsafe_log_sse2(val);
        _mm_storeu_pd(sse2_res+k, val);
        BOOST_CHECK(uni_res[k] == sse2_res[k]);
        BOOST_CHECK(uni_res[k+1] == sse2_res[k+1]);
    }
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    double avx2_res[8];
    for(unsigned int k = 0; k < 8; k+=4){
        __m256d val = _mm256_loadu_pd(values + k);
        val = unsafe_log_avx2(val);
        _mm256_storeu_pd(avx2_res+k, val);
    }
    for(unsigned int k = 0; k < 8; k++){
        BOOST_CHECK(uni_res[k] == avx2_res[k]);
    }
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    double avx512_res[8];
    __m512d val512 = _mm512_loadu_pd(values);
    val512 = unsafe_log_avx512(val512);
    _mm512_storeu_pd(avx512_res, val512);
    for(unsigned int k = 0; k < 8; k++){
        BOOST_CHECK(uni_res[k] == avx512_res[k]);
    }
#endif
}