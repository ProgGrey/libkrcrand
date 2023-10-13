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
    Xoshiro256mmUniversal gen_u[4];
    Xoshiro256mmSSE2 gen_sse[2];
    Xoshiro256mmAVX2 gen_avx;
    Xoshiro256mmState state_u = st;
    for(unsigned int k = 0; k < 4; k++){
        state_u = gen_u[k].set_state(state_u);
    }
    state_u = st;
    for(unsigned int k = 0; k < 2; k++){
        state_u = gen_sse[k].set_state(state_u);
    }
    gen_avx.set_state(st);
    unsigned int u_buf[4];
    unsigned int sse2_buf[4];
    unsigned int avx2_buf[4];
    for(unsigned int k = 0; k < Generator_Buff_Size*8; k++){
        for(unsigned int j =0; j < 4; j++){
            u_buf[j] = gen_u[j]();
        }
        for(unsigned int j = 0; j < 2; j++){
            for(unsigned int i = 0; i < 2; i++){
                sse2_buf[j*2+i] = gen_sse[j]();
            }
        }
        for(unsigned int j =0; j < 4; j++){
            avx2_buf[j] = gen_avx();
        }
        for(unsigned int j = 0; j < 4; j++){
            BOOST_CHECK(u_buf[j] == sse2_buf[j]);
        }
        for(unsigned int j = 0; j < 4; j++){
            BOOST_CHECK(u_buf[j] == avx2_buf[j]);
        }
    }
    cout << endl;
}