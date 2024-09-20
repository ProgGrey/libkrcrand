#pragma once

#include "linearlyApproximatedDistribution.hpp"
#include "exponentialDistribution.hpp"

namespace krcrand{

template<typename GenType>
class SplitedDistribution:public LAD<GenType, true, true>{
    
private:
    double Minv;
    double q_gap;
protected:
    double lambda, s, p;
    bool use_base_dist;

    ExponentialDistribution<GenType> exp_phase;
    GenType exp_prob;


    virtual double left_aprox_qf(double u) override final
    {
        return qf(u);
    }
    virtual double right_aprox_qf(double u) override final
    {
        return qf(u);
    }

    virtual double left_max(double x) override final
    {
        return 1;
    }

    virtual double right_max(double x) override final
    {
        return 1;
    }

    virtual double left_min(double x) override final
    {
        return 1;
    }
    virtual double right_min(double x) override final
    {
        return 1;
    }

    virtual double left_aprox_pdf(double x)
    {
        return pdf(x);
    }
    virtual double right_aprox_pdf(double x)
    {
        return pdf(x);
    }

    virtual double base_pdf(double x) = 0;
    virtual double base_cdf(double x) = 0;
    virtual double base_qf(double q) = 0;
   
    virtual double base_gen() = 0;
    virtual GenType::GeneratorStateType init_base(GenType::GeneratorStateType gs) = 0;

    virtual double computate_p(void) = 0;
    virtual double computate_s(void) = 0;
    virtual double computate_lambda(void) = 0;
   
    GenType::GeneratorStateType init_splited(GenType::GeneratorStateType gs)
    {
        this->p = computate_p();
        if(p > 0.0){
            gs = exp_prob.set_state(gs);
            this->lambda = computate_lambda();
            exp_phase = ExponentialDistribution<GenType>(lambda, gs);
            this->s = computate_s();
            q_gap = cdf(s);
            if (this->p >= 0.2){
                // use LAD aproximation. 
                gs = this->init_lad(gs);
                use_base_dist = false;
            } else{
                // use generator from base distribution
                Minv = 1.0-p;
                gs = this->u_generator.set_state(gs);
                gs = init_base(gs);
                use_base_dist = true;
            }
        } else{
            use_base_dist = true;
        }
        return gs;
    }
public:
    virtual double pdf(double x) override final
    {
        if(x <=s ){
            return base_pdf(x)/(1.0-p);
        } else{
            return (base_pdf(x) - p*lambda*exp(-lambda*(x-s)))/(1.0-p);
        }
    }

    double cdf(double x)
    {
        if(x <=s ){
            return base_cdf(x)/(1.0-p);
        } else{
            return (base_cdf(x) - p*(1-exp(-lambda*exp(x-s))))/(1.0-p);
        }
    }

    virtual double qf(double q) override final
    {
        if(q <= q_gap){
            return base_qf(q*(1.0-p));
        } else{
            double x = base_qf(q);
            double prev_x;
            do{
                prev_x = x;
                x = x - (cdf(x)-q)/pdf(x);
            }while(abs((prev_x - x)/x) > DBL_EPSILON);
            return x;
        }
    }

    double get_s(void)
    {
        return s;
    }

    double get_p(void)
    {
        return p;
    }

    double get_lambda(void)
    {
        return lambda;
    }

    double operator()(bool &is_exp)
    {
        if(uniform01_exclude0(exp_phase()) <= this->p){
            is_exp = true;
            return exp_phase() + s;
        }else{
            is_exp = false;
            if(use_base_dist){
                double x;
                do{
                    x = base_gen();
                }while(pdf(x)/base_pdf(x)*Minv < uniform01(this->u_generator()));
                return x;
            }else{
                return this->ladgen();
            }
        }
    }
};
}
