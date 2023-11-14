#include "platform.hpp"
#include "math.hpp"
#include <cmath>
#include <typeinfo>

namespace krcrand{

template<typename GenType> class GammaDistributionSplited
{
private:
    double x0, lambda, p;
    // Shifted exp dist:
    double exp_lam_x0, lam_m1;

    GenType generator_base;
    GenType generator_exp;
    GenType generator_gamma;
    GenType generator_p;
    double ExpGen()
    {
        return x0-lam_m1*log(exp_lam_x0 - uniform01_exclude0(generator_exp()));
    }
    // for helping dist
    double hd_add, hd_pow, hd_mult, hd_max;
    double help_dist()
    {
       return fma(hd_mult, unsafe_lambertw1(hd_pow*(-uniform01_exclude0(generator_base()))), hd_add);
    }

    double e_bm1, e_pow;
    double exp_pdf(double x)
    {
        return exp(x*e_bm1)*e_pow;
    }

    double g_mult;
    double gamma_pdf(double x)
    {
        return g_mult*pow(x,alpha - 1)*exp(x*e_bm1);
    }

    double pdf(double x)
    {
        if(x <= x0){
            return gamma_pdf(x)/(1-p);
        } else{
            return fma(exp_pdf(x), -p, gamma_pdf(x))/(1-p);
        }
    }

    double rat_e_add, rat_e_m, rat_e_den;
    double alpha_addm1, rat_num_mul;
    double rat(double x)
    {
        if(x <= x0){
            double num = exp(fma(x,rat_e_m,rat_e_add)*rat_e_den)
            num *= pow(x,alpha_addm1)*rat_num_mul;
            num/= fma(alpha_addm1, x, 1);
        } else{

        }
    }
    
    double max_value_m1;
    inline double test(double x)
    {
        return max_value_m1*rat(x);
    }

    
    void init(double alpha, double beta)
    {
        if(alpha <= 1.0 || beta <= 0.0){
            throw;
        }
        // Exp generator
        x0 = fma(alpha, beta, -beta);
        lambda = 1.0/beta;
        p = pow(alpha-1, alpha-1)/(exp(alpha-1)*tgamma(alpha));
        exp_lam_x0 = exp(lambda * x0);
        lam_m1 = 1/lambda;
        // Generator for helping dist
        double num = fma(fma(alpha, alpha, -alpha),beta,fma(beta*beta, fma(2.0, alpha, -2.0), 2.0));
        double tmp = fma(-2.0,alpha, 2.0);
        hd_add = num/tmp;
        double den = fma(2.0, beta, alpha)*fma(alpha,beta,-beta);
        hd_pow = num/den;
        hd_mult = den/tmp;
        hd_pow = hd_pow * exp(-hd_pow);
        // pdfs
        e_bm1 = -1.0/beta;
        e_pow = -e_bm1*exp(alpha - 1.0);
        g_mult = 1.0/(pow(beta,alpha)*tgamma(alpha));
        //Test function
        alpha_addm1 = alpha - 1;
        rat_e_m = -a - 2*b + 2;
        rat_e_add = (2*a - 2)*b^2 + (a^2 - a)*b;
        rat_e_den = 1/(b*(2*b + a));
        rat_num_mul = (2*b + a)*(a - 1)*b^(1 - a)*(a^2*b + (2*b^2 - b)*a - 2*b^2 + 2)/(4*(exp(a - 1)*GAMMA(a)*(a - 1) - (a - 1)^a));


    }
public:
    double operator()(){
        double x;
        do{
            x =  help_dist();
        }while(test(x) < uniform01(generator_p()));
        return x;
    }
    explicit GammaDistributionSplited(double alpha, double beta){
        pos = Generator_Buff_Size;
        generator_base.seed(0);
    }
    GammaDistributionSplited(uint64_t seed, double alpha, double beta){
        pos = Generator_Buff_Size;
        generator_base.seed(seed);
    }
};

}