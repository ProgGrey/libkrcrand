#pragma once
#include "platform.hpp"
#include "math.hpp"
#include <typeinfo>
#include <type_traits>
#include "generators.hpp"
#include <cmath>
#include <cfloat>
#include <vector>

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#endif

//#ifdef LIBKRCRAND_ENABLE_AVX2
#include <immintrin.h>
//#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
#endif

#pragma GCC diagnostic ignored "-Wignored-attributes"
namespace krcrand{

template<typename GenType, bool is_left = true, bool is_right = true> class LAD{
    private:
    double *a = nullptr, *b = nullptr, *c = nullptr;//a, b and c coefficients of 0.5*a*x^2+b*x+c.
    double *x = nullptr, *y = nullptr;
    double M = 1;// Normalising constant

    protected:
    
    GenType lad_generator;
    GenType u_generator;

    virtual double left_aprox_qf(double u) = 0;
    virtual double right_aprox_qf(double u) = 0;

    DECL_KRCRAND_ALIGN double buf_res[Generator_Buff_Size];
    DECL_KRCRAND_ALIGN double buf_u[Generator_Buff_Size];
    uint8_t pos[Generator_Buff_Size];
    unsigned int buf_pos = Generator_Buff_Size;

    // Different methods for fill result buffer
    template<typename T = GenType>
    typename std::enable_if<std::is_same<decltype(std::declval<T>().gen()), uint64_t>::value>::type fill() {
        //std::cout << "uint64_t method\n";
        for(unsigned int  k = 0; k < Generator_Buff_Size; k++){
            buf_u[k] = uniform01(u_generator());
            uint64_t u_r = lad_generator();
            double u = uniform01_exclude01(u_r);
            //return qf(u);
            // First byte equals number of interval
            uint16_t p = u_r >> tbl_shift;
            pos[k] = p;
            if(is_left && (p == 0)){
                buf_res[k] = left_aprox_qf(u);
            } else if (is_right && (p == 255)){
                buf_res[k] = right_aprox_qf(u);
            }else{
                if(is_left){
                    p--;
                }
                double cc = c[p] - u;
                //__m256d tmpv= unsafe_p2_half_solve(_mm256_set1_pd(a[p]), _mm256_set1_pd(b[p]), _mm256_set1_pd(cc), _mm256_set1_pd(x[p]), _mm256_set1_pd(x[p+1]));
                //buf_res[k] = ((double*)(&tmpv))[0];
                buf_res[k] = unsafe_p2_half_solve(a[p], b[p], cc, x[p], x[p+1]);
            }
        }
        buf_pos = 0;
    }

    #ifdef LIBKRCRAND_ENABLE_SSE2
    template<typename T = GenType>
    typename std::enable_if<std::is_same<decltype(std::declval<T>().gen()), __m128i>::value>::type fill() {
        //std::cout << "__m128i method\n";
        for(unsigned int  k = 0; k < Generator_Buff_Size; k+=2){
            //this->buf_u[k] = uniform01(u_generator());
            //this->buf_u[k+1] = uniform01(u_generator());
            _mm_store_pd(this->buf_u + k, uniform01(u_generator.gen()));
            //Generate u~U(0,1)
            __m128i u_r = lad_generator.gen();
            __m128d u = uniform01_exclude01(u_r);
            // Extract interval number
            __m128i p = _mm_srli_epi64(u_r, tbl_shift);
            //alignas(16) uint64_t p_b[2];
            //_mm_store_si128(reinterpret_cast<__m128i*>(p_b), p);
            uint64_t *p_b = reinterpret_cast<uint64_t*>(&p);
            // We can't use intrict because pos is uint8 array
            pos[k] = p_b[0];
            pos[k+1] = p_b[1];
            __m128i cond = _mm_set1_epi64x(0);
            if(is_right){
                // we can use it because p <= 255
                cond = _mm_cmpeq_epi32(p, _mm_set1_epi64x(255));
            }
            if(is_left){
                // we can use it because p <= 255
                cond = _mm_or_si128(cond, _mm_cmpeq_epi32(p, _mm_set1_epi64x(0)));
                p = _mm_sub_epi64(p, _mm_set1_epi64x(1));
            }
            // Just andnot - zero index is always valid and we ignore this computation results
            p = _mm_andnot_si128(cond, p);
            // Write coefs to registers
            __m128d av, bv, cv, x1v, x2v;
            av = _mm_set_pd(a[p_b[1]], a[p_b[0]]);
            bv = _mm_set_pd(b[p_b[1]], b[p_b[0]]);
            cv = _mm_set_pd(c[p_b[1]], c[p_b[0]]);
            x1v = _mm_set_pd(x[p_b[1]], x[p_b[0]]);
            x2v = _mm_set_pd(x[p_b[1] + 1], x[p_b[0] + 1]);
            //solve equation
            cv = _mm_sub_pd(cv, u);
            __m128d res = unsafe_p2_half_solve(av, bv, cv, x1v, x2v);
            _mm_store_pd(buf_res + k, res);
            // Fix left and right bounds
            #ifdef LIBKRCRAND_ENABLE_SSE4_1
            if(!_mm_testz_si128(cond, _mm_set1_epi64x(0xFF))){
            #else
            if(pos[k] ==0 || pos[k] == 255 || pos[k+1] == 0 || pos[k+1] == 255){
            #endif
                //cond = _mm_and_si128(cond, _mm_set1_epi64x(0xFF));
                //std::cout << (int)pos[k] << ' ' << (int)pos[k + 1] << '\n';
                alignas(16) double u_b[2];
                _mm_store_pd(u_b, u);
                for(unsigned int j = 0; j < 2; j++){
                    if(is_left && (pos[k + j] == 0)){
                        buf_res[k + j] = left_aprox_qf(u_b[j]);
                    } else if (is_right && (pos[k + j] == 255)){
                        buf_res[k + j] = right_aprox_qf(u_b[j]);
                    }
                }
            }
        }
        buf_pos = 0;
    }
    #endif

    #ifdef LIBKRCRAND_ENABLE_AVX2
    template<typename T = GenType>
    typename std::enable_if<std::is_same<decltype(std::declval<T>().gen()), __m256i>::value>::type fill() {
        for(unsigned int  k = 0; k < Generator_Buff_Size; k+=4){
            //this->buf_u[k] = uniform01(u_generator());
            //this->buf_u[k+1] = uniform01(u_generator());
            _mm256_store_pd(this->buf_u + k, uniform01(u_generator.gen()));
            //Generate u~U(0,1)
            __m256i u_r = lad_generator.gen();
            __m256d u = uniform01_exclude01(u_r);
            // Extract interval number
            __m256i p = _mm256_srli_epi64(u_r, tbl_shift);
            uint64_t *p_b = reinterpret_cast<uint64_t*>(&p);
            // We can't use intrict because pos is uint8 array
            pos[k] = p_b[0];
            pos[k+1] = p_b[1];
            pos[k+2] = p_b[2];
            pos[k+3] = p_b[3];
            __m256i cond = _mm256_set1_epi64x(0);
            if(is_right){
                cond = _mm256_cmpeq_epi64 (p, _mm256_set1_epi64x(255));
            }
            if(is_left){
                // we can use it because p <= 255
                cond = _mm256_or_si256(cond, _mm256_cmpeq_epi64(p, _mm256_set1_epi64x(0)));
                p = _mm256_sub_epi64(p, _mm256_set1_epi64x(1));
            }
            // Just andnot - zero index is always valid and we ignore this computation results
            p = _mm256_andnot_si256(cond, p);
            // Write coefs to registers
            __m256d av, bv, cv, x1v, x2v;
            av = _mm256_set_pd(a[p_b[3]], a[p_b[2]], a[p_b[1]], a[p_b[0]]);
            bv = _mm256_set_pd(b[p_b[3]], b[p_b[2]], b[p_b[1]], b[p_b[0]]);
            cv = _mm256_set_pd(c[p_b[3]], c[p_b[2]], c[p_b[1]], c[p_b[0]]);
            x1v = _mm256_set_pd(x[p_b[3]], x[p_b[2]], x[p_b[1]], x[p_b[0]]);
            x2v = _mm256_set_pd(x[p_b[3] + 1], x[p_b[2] + 1], x[p_b[1] + 1], x[p_b[0] + 1]);
            //solve equation
            cv = _mm256_sub_pd(cv, u);
            __m256d res = unsafe_p2_half_solve(av, bv, cv, x1v, x2v);
            _mm256_store_pd(buf_res + k, res);
            // Fix left and right bounds
            if(!_mm256_testz_si256(cond, _mm256_set1_epi64x(0xFF))){
                //std::cout << (int)pos[k] << ' ' << (int)pos[k + 1] << '\n';
                alignas(32) double u_b[4];
                _mm256_store_pd(u_b, u);
                for(unsigned int j = 0; j < 4; j++){
                    if(is_left && (pos[k + j] == 0)){
                        buf_res[k + j] = left_aprox_qf(u_b[j]);
                    } else if (is_right && (pos[k + j] == 255)){
                        buf_res[k + j] = right_aprox_qf(u_b[j]);
                    }
                }
            }
        }
        buf_pos = 0;
    }
    #endif

    #ifdef LIBKRCRAND_ENABLE_AVX512F
    template<typename T = GenType>
    typename std::enable_if<std::is_same<decltype(std::declval<T>().gen()), __m512i>::value>::type fill() {
        for(unsigned int  k = 0; k < Generator_Buff_Size; k+=8){
            //this->buf_u[k] = uniform01(u_generator());
            //this->buf_u[k+1] = uniform01(u_generator());
            _mm512_store_pd(this->buf_u + k, uniform01(u_generator.gen()));
            //Generate u~U(0,1)
            __m512i u_r = lad_generator.gen();
            __m512d u = uniform01_exclude01(u_r);
            // Extract interval number
            __m512i p = _mm512_srli_epi64(u_r, tbl_shift);
            uint64_t *p_b = reinterpret_cast<uint64_t*>(&p);
            // Move lower 1 byte in all elements to first element in p_wb
            __m128i p_wb =  _mm512_cvtepi64_epi8(p);
            _mm_storel_epi64(reinterpret_cast<__m128i*>(pos + k), p_wb);

            uint8_t cond = 0;
            if(is_right){
                cond = _mm512_cmp_epu64_mask(p, _mm512_set1_epi64(255), _MM_CMPINT_EQ);
            }
            if(is_left){
                // we can use it because p <= 255
                cond = cond | _mm512_cmp_epu64_mask(p, _mm512_set1_epi64(0), _MM_CMPINT_EQ);
                p = _mm512_sub_epi64(p, _mm512_set1_epi64(1));
            }
            // Set out range values to zero
            p = _mm512_mask_mov_epi64(p, cond, _mm512_set1_epi64(0x0));
            // Write coefs to registers
            __m512d av, bv, cv, x1v, x2v;
            av = _mm512_set_pd(a[p_b[7]], a[p_b[6]], a[p_b[5]], a[p_b[4]], a[p_b[3]], a[p_b[2]], a[p_b[1]], a[p_b[0]]);
            bv = _mm512_set_pd(b[p_b[7]], b[p_b[6]], b[p_b[5]], b[p_b[4]], b[p_b[3]], b[p_b[2]], b[p_b[1]], b[p_b[0]]);
            cv = _mm512_set_pd(c[p_b[7]], c[p_b[6]], c[p_b[5]], c[p_b[4]], c[p_b[3]], c[p_b[2]], c[p_b[1]], c[p_b[0]]);
            x1v = _mm512_set_pd(x[p_b[7]], x[p_b[6]], x[p_b[5]], x[p_b[4]], x[p_b[3]], x[p_b[2]], x[p_b[1]], x[p_b[0]]);
            x2v = _mm512_set_pd(x[p_b[7] + 1], x[p_b[6] + 1], x[p_b[5] + 1], x[p_b[4] + 1], x[p_b[3] + 1], x[p_b[2] + 1], x[p_b[1] + 1], x[p_b[0] + 1]);
            //solve equation
            cv = _mm512_sub_pd(cv, u);
            __m512d res = unsafe_p2_half_solve(av, bv, cv, x1v, x2v);
            _mm512_store_pd(buf_res + k, res);
            // Fix left and right bounds
            if(cond){
                alignas(64) double u_b[8];
                _mm512_store_pd(u_b, u);
                for(unsigned int j = 0; j < 8; j++){
                    if(is_left && (pos[k + j] == 0)){
                        buf_res[k + j] = left_aprox_qf(u_b[j]);
                    } else if (is_right && (pos[k + j] == 255)){
                        buf_res[k + j] = right_aprox_qf(u_b[j]);
                    }
                }
            }
        }
        buf_pos = 0;
    }
    #endif


    inline double help_dist(uint8_t &pos, double &u)
    {
        if(buf_pos >= Generator_Buff_Size){
            fill();
        }
        double res = buf_res[buf_pos];
        pos = this->pos[buf_pos];
        u = this->buf_u[buf_pos];
        buf_pos++;
        return res;
    }

    virtual double left_aprox_pdf(double x) = 0;
    virtual double right_aprox_pdf(double x) = 0;

    double pdf_approx(double x, uint8_t pos)
    {
        if(is_left && (pos == 0)){
            return left_aprox_pdf(x);
        } else if (is_right && (pos == 255)){
            return right_aprox_pdf(x);
        } else {
            uint8_t p = (is_left ? pos - 1 : pos);
            return a[p]*x + b[p];
        }
    }

    // Points where maximum of f(x)/g(x) is possible.
    virtual double left_max(double x) = 0;
    virtual double med_maxmin(double x1, double x2, double a, double b, double &y) = 0;
    virtual double right_max(double x) = 0;
    virtual double left_min(double x) = 0;
    virtual double right_min(double x) = 0;

    GenType::GeneratorStateType init_gens(GenType::GeneratorStateType gs)
    {
        gs = lad_generator.set_state(gs);
        return u_generator.set_state(gs);
    }

    GenType::GeneratorStateType init_lad(GenType::GeneratorStateType gs)
    {
        a = new double[tbl_size];
        b = new double[tbl_size];
        c = new double[tbl_size];
        x = new double[tbl_size + 1];
        y = new double[256];
        int shift = (is_left ? 1 : 0);
        double x2 = qf(static_cast<double>(shift)/256.0);
        if(is_left){
            y[0] = left_min(x2);
        }
        for(unsigned  int k = 0; k < tbl_size; k++){
            double x1 = x2;
            x2 = qf(static_cast<double>(k + shift + 1)/256.0);
            double f1 = pdf(x1);
            double f2 = pdf(x2);
             //*/
            /* Linear aproximation g(x) of PDF f(x) and quadratic aproximation G(x) of CDF F(X) with folowing properties:
            1) int_a^b f(x) = int_a^b g(x)
            2) f(a)/g(a) = f(b)/g(b)
            3) F(a) = G(a)
            4) F(b) = G(b)
            */
            double t = 2.0/256.0/((f1+f2)*(x1-x2)*(x1-x2));
            a[k] = (f2-f1)*t;
            b[k] = (x2*f1-x1*f2)*t;
            //*/
            /* Linear aproximation g(x) of PDF f(x) and quadratic aproximation G(x) of CDF F(X) with folowing properties:
            1) int_a^b f(x) = int_a^b g(x)
            2) max|f(x) - g(x)| -> min for all x in [a,b] with property 1)
            3) F(a) = G(a)
            4) F(b) = G(b)
            The following solution for a and b was obtained by condition like Chebyshev alterance theorem.
            //*/
            /*
            a[k] = (f2-f1)/(x2-x1);
            b[k] = 1.0/256.0/(x2-x1) - a[k]*(x1+x2)*0.5;
            // Fix negative density
            if(a[k]*x1 + b[k] <= 0.0){
                double t1 = 2.0/256.0/((x1-x2)*(x1-x2));
                double t2 = f1/(x2-x1);
                a[k] = t1-2*t2;
                b[k] = (x1+x2)*t2-x1*t1;
            } else if(a[k]*x2 + b[k] <= 0.0){
                double t1 = 2.0/256.0/((x1-x2)*(x1-x2));
                double t2 = f2/(x2-x1);
                a[k] = 2*t2-t1;
                b[k] = x2*t1 - (x1+x2)*t2;
            }
            //*/
            c[k] = static_cast<double>(k + shift)/256.0 - 0.5*a[k]*x1*x1 - b[k]*x1;
            x[k] = x1;
            //y[k + shift] = pdf_min(x1, x2, k + shift);
            //M = std::max(M, f1/(a[k]*x1+b[k]));
            double yy1 = f2/(a[k]*x2+b[k]);
            double yy2;
            double yy3 = med_maxmin(x1, x2, a[k], b[k], yy2);
            M = std::max(M, std::max(yy1, yy3));
            y[k + shift] = std::min(yy1, yy2);
        }
        x[tbl_size] = x2;
        if(is_left){
            M = std::max(M, left_max(x[0]));
        }
        if(is_right){
            y[255] = right_min(x2);
            M = std::max(M, right_max(x[0]));
        }
        for(unsigned int  k =0; k <= 255; k++){
            y[k] /=M;
        }
        /*
        std::cout << 'M' << M << '\n';
        std::cout.precision(16);
        std::cout << "a = c(";
        for(int  k = 0; k < tbl_size; k++){
            std::cout << a[k] << ',';
            if(k % 10 == 0){
                std::cout << '\n';
            }
        }
        std::cout << "\b)\n b = c(";
        for(int  k = 0; k < tbl_size; k++){
            std::cout << b[k] << ',';
            if(k % 10 == 0){
                std::cout << '\n';
            }
        }
        std::cout << "\b)\n cc = c(";
        for(int  k = 0; k < tbl_size; k++){
            std::cout << c[k] << ',';
            if(k % 10 == 0){
                std::cout << '\n';
            }
        }
        std::cout << "\b)\n";//*/
        return init_gens(gs);
    }

    public:
    const unsigned int tbl_size = 256 - (is_left ? 1 : 0) - (is_right ? 1 : 0);
    const unsigned int tbl_shift = 56;//64-log2(tbl_size)

    virtual double pdf(double x) = 0;
    virtual double pdf_fast(double x) = 0;
    virtual double qf(double x) = 0;

    double ladgen(){
        //return qf(uniform01_exclude01(lad_generator()));
        double x,u, tmp;
        uint8_t pos;
        do{
            x = help_dist(pos, u);
            //u = uniform01(u_generator());
            tmp = M*pdf_approx(x, pos);
        } while((y[pos] < u) && (pdf_fast(x)/(tmp) < u));
        return x;
        //return(qf(uniform01_exclude01(lad_generator())));
    }

    ~LAD()
    {
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] x;
        delete [] y;
        /*
        a = nullptr;
        b = nullptr;
        c = nullptr;
        x = nullptr;
        //*/
    }
};

    
}