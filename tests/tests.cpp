#include "../headers/libkrcrand.hpp"
#include "../headers/math.hpp"
#include "../headers/exponentialDistribution.hpp"

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cmath>

using namespace std;
using namespace krcrand;

BOOST_AUTO_TEST_CASE(constant_tests)
{
    BOOST_CHECK(sizeof(uint64_t) == sizeof(double));
    BOOST_CHECK(Generator_Buff_Size*sizeof(uint64_t) % sizeof(uint64_t) == 0);
#ifdef LIBKRCRAND_ENABLE_SSE2
    BOOST_CHECK(Generator_Buff_Size*sizeof(uint64_t) % sizeof(__m128i) == 0);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    BOOST_CHECK(Generator_Buff_Size*sizeof(uint64_t) % sizeof(__m256i) == 0);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    BOOST_CHECK(Generator_Buff_Size*sizeof(uint64_t) % sizeof(__m512i) == 0);
#endif
}


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
    for(unsigned int k = 0; k < 100*Generator_Buff_Size; k++){
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
}


BOOST_AUTO_TEST_CASE(math_tests)
{
    const double values[8] = {0.1, 2, 3.5, 1, 0.000001, 12, 7643921, 2.2250738585072009e-318};
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
        val = unsafe_log(val);
        _mm_storeu_pd(sse2_res+k, val);
        BOOST_CHECK(uni_res[k] == sse2_res[k]);
        BOOST_CHECK(uni_res[k+1] == sse2_res[k+1]);
    }
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    double avx2_res[8];
    for(unsigned int k = 0; k < 8; k+=4){
        __m256d val = _mm256_loadu_pd(values + k);
        val = unsafe_log(val);
        _mm256_storeu_pd(avx2_res+k, val);
    }
    for(unsigned int k = 0; k < 8; k++){
        BOOST_CHECK(uni_res[k] == avx2_res[k]);
    }
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    double avx512_res[8];
    __m512d val512 = _mm512_loadu_pd(values);
    val512 = unsafe_log(val512);
    _mm512_storeu_pd(avx512_res, val512);
    for(unsigned int k = 0; k < 8; k++){
        BOOST_CHECK(uni_res[k] == avx512_res[k]);
    }
#endif
}

BOOST_AUTO_TEST_CASE(exponential_distribution_tests)
{
    Xoshiro256mmState st;
    st.seed(1);
    ExponentialDistribution<Xoshiro256mmUniversalStable, 0> exp_u(9.8);
    exp_u.generator.set_state(st);
    ExponentialDistribution<Xoshiro256mmUniversalStable, 1> exp_u_n(9.8);
    exp_u_n.generator.set_state(st);
#ifdef LIBKRCRAND_ENABLE_SSE2
    ExponentialDistribution<Xoshiro256mmSSE2stable, 0> exp_sse2(9.8);
    exp_sse2.generator.set_state(st);
    ExponentialDistribution<Xoshiro256mmSSE2stable, 1> exp_sse2_n(9.8);
    exp_sse2_n.generator.set_state(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    ExponentialDistribution<Xoshiro256mmAVX2stable, 0> exp_avx2(9.8);
    exp_avx2.generator.set_state(st);
    ExponentialDistribution<Xoshiro256mmAVX2stable, 1> exp_avx2_n(9.8);
    exp_avx2_n.generator.set_state(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    ExponentialDistribution<Xoshiro256mmAVX512Fstable, 0> exp_avx512(9.8);
    exp_avx512.generator.set_state(st);
    ExponentialDistribution<Xoshiro256mmAVX512Fstable, 1> exp_avx512_n(9.8);
    exp_avx512_n.generator.set_state(st);
#endif
    unsigned int N = 10;
    double mx = 0,my = 0,mxy = 0;
    for(unsigned int k = 0; k < N; k++){
        double val = exp_u();
        mx += val/N;
        double val_n = exp_u_n();
        my += val_n/N;
        mxy += val*val_n/N;
        //cout << val << '|' << val_n << ',';
#ifdef LIBKRCRAND_ENABLE_SSE2
        BOOST_CHECK(fabs(val - exp_sse2())/val < 1e-15);
        BOOST_CHECK(fabs(val_n - exp_sse2_n())/val_n  < 1e-15);
        //cout << exp_sse2_n() << ',';
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
        BOOST_CHECK(fabs(val - exp_avx2())/val < 1e-15);
        BOOST_CHECK(fabs(val_n - exp_avx2_n())/val_n < 1e-15);
        //cout <<  exp_avx2_n() << '\n';
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
        BOOST_CHECK(val == exp_avx512());
        BOOST_CHECK(val_n == exp_avx512_n());
        //cout <<  exp_avx512_n() << endl;
#endif
    }
    // Antithetic variates test:
    BOOST_CHECK((mxy-mx*my) < -0.006);
}
//*/