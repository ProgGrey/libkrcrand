#include "../headers/libkrcrand.hpp"
#include "../headers/math.hpp"
#include "../headers/exponentialDistribution.hpp"
#include "../headers/gammaDistribution.hpp"

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

BOOST_AUTO_TEST_CASE(log_tests)
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

BOOST_AUTO_TEST_CASE(p2_solve_tests)
{
    DECL_KRCRAND_ALIGN double a[] = { 2,   -0,  2,   12, 4,  2, 2, 12};
    DECL_KRCRAND_ALIGN double b[] = { 1,   1,  0,   1,  0,  1, 0, 1};
    DECL_KRCRAND_ALIGN double c[] = {-1,   2, -1,   0,  0, -1,-1, 0};
    DECL_KRCRAND_ALIGN double x1[] = {0, -2.5, 0, -0.2, -1,-2,-2,-0.1};
    DECL_KRCRAND_ALIGN double x2[] = {1, -1.5, 2, -0.1,  1,-1, 0, 0.1};
    double res[8];
    for(unsigned int k = 0; k < 8; k++){
        res[k] = unsafe_p2_half_solve(a[k], b[k], c[k], x1[k], x2[k]);
    }
    #ifdef LIBKRCRAND_ENABLE_SSE2
    alignas(16) double res_sse[2];
    for(unsigned int  i = 0; i < 8; i++){
        for(unsigned int  j = 0; j < 8; j++){
            __m128d a_sse =  _mm_set_pd(a[j], a[i]);
            __m128d b_sse =  _mm_set_pd(b[j], b[i]);
            __m128d c_sse =  _mm_set_pd(c[j], c[i]);
            __m128d x1_sse =  _mm_set_pd(x1[j], x1[i]);
            __m128d x2_sse =  _mm_set_pd(x2[j], x2[i]);
            __m128d r_sse = unsafe_p2_half_solve(a_sse, b_sse, c_sse, x1_sse, x2_sse);
            _mm_store_pd(res_sse, r_sse);
            //cout << res[j] << ' ' << res[i] << '\n';
            BOOST_CHECK(res[i] == res_sse[0]);
            BOOST_CHECK(res[j] == res_sse[1]);    
        }
    }
    #endif
    #ifdef LIBKRCRAND_ENABLE_AVX2
    alignas(32) double res_avx[4];
    for(unsigned int  l = 0; l < 8; l++){
        for(unsigned int  k = 0; k < 8; k++){
            for(unsigned int  i = 0; i < 8; i++){
                for(unsigned int  j = 0; j < 8; j++){
                    __m256d a_avx =  _mm256_set_pd(a[j], a[i], a[k], a[l]);
                    __m256d b_avx =  _mm256_set_pd(b[j], b[i], b[k], b[l]);
                    __m256d c_avx =  _mm256_set_pd(c[j], c[i], c[k], c[l]);
                    __m256d x1_avx =  _mm256_set_pd(x1[j], x1[i], x1[k], x1[l]);
                    __m256d x2_avx =  _mm256_set_pd(x2[j], x2[i], x2[k], x2[l]);
                    __m256d r_avx = unsafe_p2_half_solve(a_avx, b_avx, c_avx, x1_avx, x2_avx);
                    _mm256_store_pd(res_avx, r_avx);
                    //cout << res[j] << ' ' << res[i] << '\n';
                    BOOST_CHECK(res[l] == res_avx[0]);
                    BOOST_CHECK(res[k] == res_avx[1]);
                    BOOST_CHECK(res[i] == res_avx[2]);
                    BOOST_CHECK(res[j] == res_avx[3]);
                }
            }
        }
    }
    #endif
    #ifdef LIBKRCRAND_ENABLE_AVX512F
    alignas(64) double res_avx512[8];
    __m512d a_avx512 =  _mm512_load_pd(a);
    __m512d b_avx512 =  _mm512_load_pd(b);
    __m512d c_avx512 =  _mm512_load_pd(c);
    __m512d x1_avx512 = _mm512_load_pd(x1);
    __m512d x2_avx512 = _mm512_load_pd(x2);
    __m512d r_avx512 = unsafe_p2_half_solve(a_avx512, b_avx512, c_avx512, x1_avx512, x2_avx512);
    _mm512_store_pd(res_avx512, r_avx512);
    for( unsigned int k =0; k < 8; k++){
        BOOST_CHECK(res[k] == res_avx512[k]);
    }
    for(unsigned int  l = 0; l < 8; l++){
        for(unsigned int  k = 0; k < 8; k++){
            for(unsigned int  i = 0; i < 8; i++){
                for(unsigned int  j = 0; j < 8; j++){
                    __m512d a_avx512 =  _mm512_set_pd(a[j], a[i], a[k], a[l], a[l], a[k], a[i], a[j]);
                    __m512d b_avx512 =  _mm512_set_pd(b[j], b[i], b[k], b[l], b[l], b[k], b[i], b[j]);
                    __m512d c_avx512 =  _mm512_set_pd(c[j], c[i], c[k], c[l], c[l], c[k], c[i], c[j]);
                    __m512d x1_avx512 =  _mm512_set_pd(x1[j], x1[i], x1[k], x1[l], x1[l], x1[k], x1[i], x1[j]);
                    __m512d x2_avx512 =  _mm512_set_pd(x2[j], x2[i], x2[k], x2[l], x2[l], x2[k], x2[i], x2[j]);
                    __m512d r_avx512 = unsafe_p2_half_solve(a_avx512, b_avx512, c_avx512, x1_avx512, x2_avx512);
                    _mm512_store_pd(res_avx512, r_avx512);
                    //cout << res[j] << ' ' << res[i] << '\n';
                    BOOST_CHECK(res[j] == res_avx512[0]);
                    BOOST_CHECK(res[i] == res_avx512[1]);
                    BOOST_CHECK(res[k] == res_avx512[2]);
                    BOOST_CHECK(res[l] == res_avx512[3]);
                    BOOST_CHECK(res[l] == res_avx512[4]);
                    BOOST_CHECK(res[k] == res_avx512[5]);
                    BOOST_CHECK(res[i] == res_avx512[6]);
                    BOOST_CHECK(res[j] == res_avx512[7]);
                }
            }
        }
    }
    #endif
}


BOOST_AUTO_TEST_CASE(exponential_distribution_tests)
{
    Xoshiro256mmState st;
    st.seed(1);
    ExponentialDistribution<Xoshiro256mmUniversalStable, 0> exp_u(9.8);
    exp_u.set_state(st);
    ExponentialDistribution<Xoshiro256mmUniversalStable, 1> exp_u_n(9.8);
    exp_u_n.set_state(st);
#ifdef LIBKRCRAND_ENABLE_SSE2
    ExponentialDistribution<Xoshiro256mmSSE2stable, 0> exp_sse2(9.8);
    exp_sse2.set_state(st);
    ExponentialDistribution<Xoshiro256mmSSE2stable, 1> exp_sse2_n(9.8);
    exp_sse2_n.set_state(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    ExponentialDistribution<Xoshiro256mmAVX2stable, 0> exp_avx2(9.8);
    exp_avx2.set_state(st);
    ExponentialDistribution<Xoshiro256mmAVX2stable, 1> exp_avx2_n(9.8);
    exp_avx2_n.set_state(st);
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    ExponentialDistribution<Xoshiro256mmAVX512Fstable, 0> exp_avx512(9.8);
    exp_avx512.set_state(st);
    ExponentialDistribution<Xoshiro256mmAVX512Fstable, 1> exp_avx512_n(9.8);
    exp_avx512_n.set_state(st);
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


// @return 0 if test is ok, 1 if KS test fails, 2 or 3 if independency tests fails
int check_distribution(const char *dist, const vector<double> &params, const vector<double> &data)
{
    size_t d_len = strlen(dist);
    const char *script = "./tests/dists.R";
    size_t s_len = strlen(script);
    char *exec_str = new char[s_len + 1 + d_len + 1];
    char *p = exec_str;
    // Copy script path
    memcpy(p, script, s_len);
    p+=s_len;
    *p= ' ';
    p++;
    // Copy dist name
    memcpy(p, dist, d_len);
    p += d_len;
    *p = '\0';
    // Write parameters
    FILE* fd = popen(exec_str, "w");
    for(double val : params){
        fprintf(fd, "%.16e ", val);
    }
    // Write data
    for(double val : data){
        fprintf(fd, "%.16e ", val);
    }
    return(pclose(fd));
}


int gamma_goodness(double alpha, double beta)
{
    Xoshiro256mm::GeneratorStateType inter(1);
    //GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, inter);
    GammaDistribution<Xoshiro256mm> gen(alpha, beta, inter);
    vector<double> params = {alpha, beta};
    vector<double> data;
    for(int k = 0; k < 1000; k++){
        data.push_back(gen());
        //std::cout << data.back() << ", ";
    }
    return check_distribution("gamma", params, data);
}


int exponential_goodness(double lambda)
{
    Xoshiro256mm::GeneratorStateType inter(1);
    //GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, inter);
    ExponentialDistribution<Xoshiro256mm> gen(lambda, inter);
    vector<double> params = {lambda};
    vector<double> data;
    for(int k = 0; k < 1000; k++){
        data.push_back(gen());
        //std::cout << data.back() << ", ";
    }
    return check_distribution("exponential", params, data);
}


BOOST_AUTO_TEST_CASE(goodness_of_fit_tests)
{
    BOOST_CHECK(exponential_goodness(0.005)== 0);
    BOOST_CHECK(exponential_goodness(0.5)== 0);
    BOOST_CHECK(exponential_goodness(1)== 0);
    BOOST_CHECK(exponential_goodness(12)== 0);
    BOOST_CHECK(exponential_goodness(10000)== 0);


    BOOST_CHECK(gamma_goodness(0.005, 10)== 0);
    BOOST_CHECK(gamma_goodness(0.5, 0.1)== 0);
    BOOST_CHECK(gamma_goodness(1, 1)== 0);
    BOOST_CHECK(gamma_goodness(2, 0.1)== 0);
    BOOST_CHECK(gamma_goodness(130, 300)== 0);
    BOOST_CHECK(gamma_goodness(1300, 0.5)== 0);

    
}