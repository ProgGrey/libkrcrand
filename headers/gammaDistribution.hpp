#include "platform.hpp"
#include "math.hpp"
#include <cmath>
#include <typeinfo>
#include "exponentialDistribution.hpp"

namespace krcrand{

template<typename GenType, typename GeneratorStateType = GenType::GeneratorStateType> class GammaDistributionSplited
{
private:
    double x0, lambda, p;
    // Shifted exp dist:
    double exp_lam_x0, lam_m1;

    GenType generator_base;
    GenType generator_exp;
    ExponentialDistribution<GenType, 0> exp_dist;
    GenType generator_gamma;
    GenType generator_p;
    double exp_gen()
    {
        return exp_dist() + x0;
        return x0-lam_m1*log(uniform01_exclude0(generator_exp()));
        return x0-lam_m1*log(exp_lam_x0 - uniform01_exclude0(generator_exp()));
    }
    // for helping dist
    double hd_mult, hd_lam_mult, hd_add;
    double help_dist()
    {
        double um1 = -uniform01_exclude0(generator_gamma());
        //std::cout << "um1 " << um1 << std::endl;
        //std::cout << "hd_lam_mult " << hd_lam_mult << std::endl;
        //std::cout << "hd_add " << hd_add << std::endl;
        //std::cout << "hd_mult " << hd_mult << std::endl;
        //hd_lam_mult+=0.00001;
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


    GeneratorStateType init_gens(GeneratorStateType state, bool is_exp)
    {
        state = generator_base.set_state(state);
        state = generator_gamma.set_state(state);
        state = generator_p.set_state(state);
        if(is_exp){
            state = generator_exp.set_state(state);
        }
        return state;
    }
    
    GeneratorStateType init(double alpha, double beta, GeneratorStateType gs)
    {
        if(alpha <= 1.0 || beta <= 0.0){
            throw;
        }
        // Exp generator
        x0 = fma(alpha, beta, -beta);
        lambda = 1.0/beta;
        gs = init_gens(gs, false);
        exp_dist = ExponentialDistribution<GenType, 0>(lambda, gs);
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
        //std::cout << "tmp " << tmp << std::endl;
        hd_lam_mult = tmp*exp(-tmp);
        if(hd_lam_mult == 0){
            throw;
        }
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
        //std::cout << max_value_m1 << std::endl;
        //std::cout << x0 << std::endl;
        max_value_m1 = 1.0/max_value_m1;
        return gs;
    }

    bool exp_flag;
public:
    uint64_t max_fails=0;
    uint64_t sum_fails=0;
    double operator()(){
        //return help_dist();
        //return exp_gen();
        double x;
        if(uniform01(generator_p()) < p){
            exp_flag = true;
            return exp_gen();
        }
        exp_flag = false;
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

    double get_x0()
    {
        return x0;
    }

    bool is_exp_phase()
    {
        return exp_flag;
    }

    GammaDistributionSplited(){
        GeneratorStateType state(0);
        init(2, 1, state);
    }

    explicit GammaDistributionSplited(double alpha, double beta){
        GeneratorStateType state(0);
        init(alpha, beta, state);
    }
    explicit GammaDistributionSplited(double alpha, double beta, uint64_t seed){
        GeneratorStateType state(seed);
        init(alpha, beta, state);
    }
    explicit GammaDistributionSplited(double alpha, double beta, GeneratorStateType state){
        init(alpha, beta, state);
    }

    GeneratorStateType set_state(GeneratorStateType state)
    {
        return init_gens(state, true);
    }

};

}