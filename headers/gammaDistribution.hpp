#include "platform.hpp"
#include "math.hpp"
#include <cmath>
#include <typeinfo>

#include <iostream>

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
    double exp_gen()
    {
        return x0-lam_m1*log(uniform01_exclude0(generator_exp()));
        return x0-lam_m1*log(exp_lam_x0 - uniform01_exclude0(generator_exp()));
    }
    // for helping dist
    double hd_mult, hd_lam_mult, hd_add;
    double help_dist()
    {
        double um1 = -uniform01_exclude0(generator_gamma());
        //std::cout << "um1 " << um1 << std::endl;
        return hd_mult*unsafe_lambertw1(um1*hd_lam_mult) + hd_add;
    }

    double left_exp_m, left_exp_a;
    double left_mull;
    double a_subs_2;
    double a_subs_1, right_x_mull, right_x_add, right_mull;
    double rat_left(double x)
    {
        double e = exp(fma(left_exp_m,x,left_exp_a));
        return left_mull*e*pow(x,a_subs_1)/(x+1.0);
    }
    double rat_right(double x)
    {
        double e = exp(fma(left_exp_m,x,left_exp_a));
        return e*right_mull*fma(pow(x, a_subs_1),right_x_mull, right_x_add)/(x+1.0);
    }
    double rat(double x)
    {
        if(x <= x0){
            return rat_left(x);
        } else{
            return rat_right(x);
        }
    }
    
    double max_value_m1;
    inline double test(double x)
    {
        return max_value_m1*rat(x);
    }

    
    void init(double alpha, double beta, uint64_t seed)
    {
        auto gs = generator_base.seed(seed);
        gs = generator_exp.set_state(gs);
        gs = generator_gamma.set_state(gs);
        generator_p.set_state(gs);
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
        double sqsqa = sqrt(sqrt(alpha));
        double sqa = sqrt(alpha);
        double a34 = sqa*sqsqa;
        hd_mult = -beta*a34;
        hd_add = hd_mult - 1.0;
        double tmp = 1.0+1.0/(a34*beta);
        hd_lam_mult = tmp*exp(-tmp);
        // Test function x<=x0:
        left_exp_m = (1.0-a34)/(a34*beta);
        left_exp_a = alpha - 1;
        a_subs_2 = alpha - 2;
        left_mull = sqa*((sqsqa*alpha-sqsqa)*pow(beta,1.0-alpha)+alpha*pow(beta,2.0-alpha)*(alpha-1.0))/
                    (exp(alpha - 1.0)*tgamma(alpha)*(alpha-1.0)-pow(alpha - 1.0, alpha));
        //Test function x > x0
        right_mull = sqa*fma(alpha,beta, sqsqa)/(exp(alpha - 1.0)*tgamma(alpha)-pow(alpha - 1.0, alpha - 1.0));
        a_subs_1 = alpha - 1.0;
        right_x_mull = pow(beta, 1.0-alpha);
        right_x_add = -pow(alpha - 1.0, alpha - 1.0);
        // Find maximum
        max_value_m1 = rat_left(0);
        tmp = rat_left(x0);
        if(tmp > max_value_m1){
            max_value_m1 = tmp;
        }
        tmp = (1.0+sqa*sqsqa*(alpha*beta-(2.0*beta+1.0)) + sqrt(1 + ((beta*(alpha*beta - 4.0*beta + 2.0)*alpha - 2.0*sqsqa*beta + (4.0*beta)*beta + 1.0)*alpha - 2.0*sqsqa)*sqa))/(2*(sqa*sqsqa-1.0));
        if(tmp > 0 && tmp < x0){
            tmp = rat_left(tmp);
            if(tmp > max_value_m1){
                max_value_m1 = tmp;
            }
        }
        std::cout << max_value_m1 << std::endl;
        std::cout << x0 << std::endl;
        max_value_m1 = 1.0/max_value_m1;
    }
public:
    uint64_t max_fails=0;
    uint64_t sum_fails=0;
    double operator()(){
        //return help_dist();
        //return exp_gen();
        double x;
        if(uniform01(generator_p()) < p){
            return exp_gen();
        }
        int fails = -1;
        do{
            x =  help_dist();
            fails++;
           //std::cout << "hd: " <<  x << std::endl;
        }while(test(x) < uniform01(generator_base()));
        //std::cout << "x: " <<  x << std::endl;
        sum_fails += fails;
        if(max_fails < fails){
            max_fails = fails;
        }
        return x;
    }
    explicit GammaDistributionSplited(double alpha, double beta){
        init(alpha, beta, 0);
    }
    GammaDistributionSplited(uint64_t seed, double alpha, double beta){
        init(alpha, beta, seed);
    }
};

}