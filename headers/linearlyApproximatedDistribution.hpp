#pragma once
#include "platform.hpp"
#include "math.hpp"
#include <typeinfo>
#include "generators.hpp"
#include <cmath>
#include <cfloat>
#include <vector>

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

    double help_dist(uint8_t &pos)
    {
        uint64_t u_r = lad_generator();
        double u = uniform01_exclude01(u_r);
        //return qf(u);
        // First byte equals number of interval
        uint16_t p = u_r >> tbl_shift;
        pos = p;
        if(is_left && (p == 0)){
            return left_aprox_qf(u);
        } else if (is_right && (p == 255)){
            return right_aprox_qf(u);
        }else{
            if(is_left){
                p--;
            }
            double cc = c[p] - u;
            /*
            if(abs(a[p]) < DBL_EPSILON){
                // Linear equation:
                return -cc/b[p];
            }else if (abs(b[p]) < DBL_EPSILON){
                // Symetric quaratic equation
                return sqrt(-2.0*cc/a[p]);
            } else if(abs(c[p]) < DBL_EPSILON){
                // Incomplete quadratic equation
                return -2.0*b[p]/a[p];
            }
            double sd = sqrt(b[p]*b[p] - 2 * a[p]*cc);
            // This root is more probable due to properties of used approximation:
            double xx = (-b[p]+sd)/a[p];
            if((x[p] <= xx) && (xx <= x[p+1])){
                return xx;
            } else{
                return (-b[p]-sd)/a[p];
            }
            //*/
            return unsafe_p2_half_solve(a[p], b[p], cc, x[p], x[p+1]);
        }
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
    virtual double med_max(double x1, double x2, double a, double b) = 0;
    virtual double right_max(double x) = 0;
    virtual double pdf_min(double x1, double x2, uint8_t p) = 0;

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
            y[0] = pdf_min(qf(0.0), x2, 0);
        }
        for(unsigned  int k = 0; k < tbl_size; k++){
            /* Linear aproximation g(x) of PDF f(x) and quadratic aproximation G(x) of CDF F(X) with folowing properties:
            1) int_a^b f(x) = int_a^b g(x)
            2) max|f(x) - g(x)| -> min for all x in [a,b] with property 1)
            3) F(a) = G(a)
            4) F(b) = G(b)
            The following solution for a and b was obtained by condition like Chebyshev alterance theorem.
            //*/
            double x1 = x2;
            x2 = qf(static_cast<double>(k + shift + 1)/256.0);
            double f1 = pdf(x1);
            double f2 = pdf(x2);
            a[k] = (f2-f1)/(x2-x1);
            b[k] = 1.0/256.0/(x2-x1) - a[k]*(x1+x2)*0.5;
            // Fix negative density
            if(a[k]*x1 + b[k] < 0.0){
                double tmp = (x1-x2);
                tmp *= tmp;
                a[k] = 2.0/256.0/tmp;
                b[k] = -2.0/256.0*x1/tmp;
            } else if(a[k]*x2 + b[k] < 0.0){
                double tmp = (x1-x2);
                tmp *= tmp;
                a[k] = -2.0/256.0/tmp;
                b[k] = 2.0/256.0*x1/tmp;
            }
            c[k] = static_cast<double>(k + shift)/256.0 - 0.5*a[k]*x1*x1 - b[k]*x1;
            x[k] = x1;
            y[k + shift] = pdf_min(x1, x2, k + shift);
            M = std::max(M, f1/(a[k]*x1+b[k]));
            M = std::max(M, f2/(a[k]*x2+b[k]));
            M = std::max(M, med_max(x1, x2, a[k], b[k]));
        }
        x[tbl_size] = x2;
        y[255] = pdf_min(x2, qf(1), 255);
        if(is_left){
            M = std::max(M, left_max(x[0]));
        }
        if(is_right){
            M = std::max(M, right_max(x[0]));
        }
        
        /*
        std::cout << 'M' << M << '\n';
        std::cout.precision(16);
        for(int  k = 0; k < tbl_size; k++){
            std::cout << a[k] << ',';
        }
        std::cout << '\n';
        for(int  k = 0; k < tbl_size; k++){
            std::cout << b[k] << ',';
        }
        std::cout << '\n';
        for(int  k = 0; k < tbl_size; k++){
            std::cout << c[k] << ',';
        }
        std::cout << '\n';//*/
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
            x = help_dist(pos);
            u = uniform01(u_generator());
            tmp = M*pdf_approx(x, pos);
        } while((y[pos]/tmp < u) && (pdf_fast(x)/(tmp) < u));
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