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
    Xoshiro256mmUniversal gen_u1;
    Xoshiro256mmUniversal gen_u2;
    Xoshiro256mmSSE2 gen_sse1;
    auto state_u = gen_u1.set_state(st);
    state_u = gen_u2.set_state(state_u);
    gen_sse1.set_state(st);
    for(unsigned int k = 0; k < Generator_Buff_Size*8; k++){
        BOOST_CHECK(gen_u1() == gen_sse1());
        BOOST_CHECK(gen_u2() == gen_sse1());
    }
    cout << endl;
}