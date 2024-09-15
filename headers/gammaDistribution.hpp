#pragma once
#include "platform.hpp"
#include "math.hpp"
#include <cmath>
#include <typeinfo>
#include "exponentialDistribution.hpp"
#include "linearlyApproximatedDistribution.hpp"

namespace krcrand{

namespace nGammaDistributionTools
{
    double computate_p(double alpha);
    double gamma_cdf(double x, double alpha, double beta);
    double gamma_qf(double p, double alpha, double beta);
    double gamma_pdf(double x, double alpha, double beta);
    double gamma_pdf_fast(double x, double alpha, double beta);
}


template<typename GenType>
class GammaDistribution final : public LAD<GenType, true, true> 
{
    using LAD<GenType, true, true>::LAD;
private:
    double alpha, beta;
    
    virtual double left_aprox_qf(double u) override
    {
        return qf(u);
    }
    virtual double right_aprox_qf(double u) override
    {
        return qf(u);
    }

    virtual double left_max(double x) override
    {
        return 1;
    }

    virtual double med_max(double x1, double x2, double a, double b) override
    {
        double xc1, xc2;
        p2_solve(a, b-a*beta*(alpha - 2), b*beta*(1.0-alpha),xc1, xc2);
        if((x1 < xc1 && xc1 < x2) && (x1 < xc2 && xc2 < x2)){
            return std::max<double>(pdf(xc1)/(a*xc1+b), pdf(xc2)/(a*xc2+b));
        } else if (x1 < xc1 && xc1 < x2){
            return pdf(xc1)/(a*xc1+b);
        } else if (x1 < xc2 && xc2 < x2){
            return pdf(xc2)/(a*xc2+b);
        } else{
            double xc = (x2+x1)*0.5;
            return pdf(xc)/(a*xc+b);
        }//*/
        return 1;
    }

    virtual double right_max(double x) override
    {
        return 1;
    }

    virtual double pdf_min(double x1, double x2, uint8_t p) override final
    {
        switch (p)
        {
        case 0:
            return 0;
            break;
        case 255:
            return 0;
            break;
        default:
            return std::min<double>(pdf(x1), pdf(x2));
            break;
        }
    }

    virtual double left_aprox_pdf(double x)
    {
        return pdf_fast(x);
    }
    virtual double right_aprox_pdf(double x)
    {
        return pdf_fast(x);
    }

    GenType::GeneratorStateType init(double alpha, double beta, GenType::GeneratorStateType gs)
    {
        this->alpha = alpha;
        this->beta = beta;
        if(alpha >= 1.0){
            return this->init_lad(gs);
        }else{
            return this->init_gens(gs);
        }
    }


public:

    virtual double pdf(double x) override
    {
        return nGammaDistributionTools::gamma_pdf(x, alpha, beta);
    }

    virtual double pdf_fast(double x) override
    {
        return nGammaDistributionTools::gamma_pdf_fast(x, alpha, beta);
    }

    virtual double qf(double x) override
    {
        return nGammaDistributionTools::gamma_qf(x, alpha, beta);
    }

    GammaDistribution(){
        using GT = GenType::GeneratorStateType;
        GT state(0);
        init(2, 1, state);
    }

    explicit GammaDistribution(double alpha, double beta){
        using GT = GenType::GeneratorStateType;
        GT state(0);
        init(alpha, beta, state);
    }

    explicit GammaDistribution(double alpha, double beta, uint64_t seed){
        using GT = GenType::GeneratorStateType;
        GT state(seed);
        init(alpha, beta, state);
    }
    explicit GammaDistribution(double alpha, double beta, GenType::GeneratorStateType &state){
        state = init(alpha, beta, state);
    }

    GenType::GeneratorStateType set_state(GenType::GeneratorStateType state)
    {
        return this->init_gens(state);
    }

    int counter = 0;
    double operator()(void)
    {
        if(alpha >= 1.0){
            return this->ladgen();
        } else{
            //GS algorithm. See J. H. Ahrens; U. Dieter. (1974). Computer methods for sampling from gamma, beta, poisson and bionomial distributions. , 12(3), 223–246. doi:10.1007/bf02293108
            double x;
            double t;
            double v;
            do{
                counter++;
                double u = uniform01_exclude0(this->u_generator());
                double b = (M_E + alpha)/M_E;
                double p = b*u;
                v = uniform01_exclude0(this->lad_generator());
                if(p > 1){
                    x = -log((b-p)/alpha);
                    t = pow(x,alpha - 1.0);
                } else {
                    x = pow(p, 1.0/alpha);
                    t = exp(-x);
                }
            }while(t < v);
            return x*beta;
            //*/
        }
    }
};

template<typename GenType>
class GammaDistributionSplited final : public LAD<GenType, true, true> 
{
    using LAD<GenType, true, true>::LAD;
private:
    double alpha, beta;
    double p, s;
    ExponentialDistribution<GenType> exp_phase;
    GenType split_gen;

    virtual double left_aprox_qf(double u) override
    {
        return qf(u);
    }
    virtual double right_aprox_qf(double u) override
    {
        return qf(u);
    }

    virtual double left_max(double x) override
    {
        return 1;
    }

    virtual double med_max(double x1, double x2, double a, double b) override
    {
        double xc1, xc2;
        p2_solve(a, b-a*beta*(alpha - 2), b*beta*(1.0-alpha),xc1, xc2);
        if((x1 < xc1 && xc1 < x2) && (x1 < xc2 && xc2 < x2)){
            return std::max<double>(pdf(xc1)/(a*xc1+b), pdf(xc2)/(a*xc2+b));
        } else if (x1 < xc1 && xc1 < x2){
            return pdf(xc1)/(a*xc1+b);
        } else if (x1 < xc2 && xc2 < x2){
            return pdf(xc2)/(a*xc2+b);
        } else{
            double xc = (x2+x1)*0.5;
            return pdf(xc)/(a*xc+b);
        }//*/
        return 1;
    }

    virtual double right_max(double x) override
    {
        return 1;
    }

    virtual double pdf_min(double x1, double x2, uint8_t p) override final
    {
        switch (p)
        {
        case 0:
            return INFINITY;
            break;
        case 255:
            return INFINITY;
            break;
        default:
            return std::min<double>(pdf(x1), pdf(x2));
            break;
        }
    }

    virtual double left_aprox_pdf(double x)
    {
        return pdf_fast(x);
    }
    virtual double right_aprox_pdf(double x)
    {
        return pdf_fast(x);
    }

    GenType::GeneratorStateType init(double alpha, double beta, GenType::GeneratorStateType gs)
    {
        this->alpha = alpha;
        this->beta = beta;
        //TODO code for lambda, p, s computations
        double lambda = 1;
        p = 0;
        s = 0;
        exp_phase = ExponentialDistribution<GenType>(lambda, gs);
        if(alpha >= 1.0){
            return this->init_lad(gs);
        }else{
            gs = this->lad_generator.set_state(gs);
            return this->u_generator.set_state(gs);
        }
    }


public:

    virtual double pdf(double x) override
    {
        return nGammaDistributionTools::gamma_pdf(x, alpha, beta);
    }

    virtual double pdf_fast(double x) override
    {
        return nGammaDistributionTools::gamma_pdf_fast(x, alpha, beta);
    }

    virtual double qf(double x) override
    {
        return nGammaDistributionTools::gamma_qf(x, alpha, beta);
    }

    GammaDistributionSplited(){
        using GT = GenType::GeneratorStateType;
        GT state(0);
        init(2, 1, state);
    }

    explicit GammaDistributionSplited(double alpha, double beta){
        using GT = GenType::GeneratorStateType;
        GT state(0);
        init(alpha, beta, state);
    }

    explicit GammaDistributionSplited(double alpha, double beta, uint64_t seed){
        using GT = GenType::GeneratorStateType;
        GT state(seed);
        init(alpha, beta, state);
    }
    explicit GammaDistributionSplited(double alpha, double beta, GenType::GeneratorStateType &state){
        state = init(alpha, beta, state);
    }

    GenType::GeneratorStateType set_state(GenType::GeneratorStateType state)
    {
        return init_gens(state);
    }

    int counter = 0;
    double operator()(void)
    {
        if(uniform01(split_gen()) < p){
            return s + exp_phase();    
        }
        if(alpha >= 1.0){
            return this->ladgen();
        } else{
            //GS algorithm. See J. H. Ahrens; U. Dieter. (1974). Computer methods for sampling from gamma, beta, poisson and bionomial distributions. , 12(3), 223–246. doi:10.1007/bf02293108
            double x;
            double t;
            double v;
            do{
                counter++;
                double u = uniform01_exclude0(this->u_generator());
                double b = (M_E + alpha)/M_E;
                double p = b*u;
                v = uniform01_exclude0(this->lad_generator());
                if(p > 1){
                    x = -log((b-p)/alpha);
                    t = pow(x,alpha - 1.0);
                } else {
                    x = pow(p, 1.0/alpha);
                    t = exp(-x);
                }
            }while(t < v);
            return x*beta;
            //*/
        }
    }
};

}