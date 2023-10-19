#include "../headers/libkrcrand.hpp"


#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace std;
using namespace krcrand;


BOOST_AUTO_TEST_CASE(Xoshiro256mm)
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