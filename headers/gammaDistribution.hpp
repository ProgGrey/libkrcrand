#pragma once
#include "platform.hpp"
#include "math.hpp"
#include <cmath>
#include <typeinfo>
#include "linearlyApproximatedDistribution.hpp"
#include "splitedDistribution.hpp"

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

    virtual double med_maxmin(double x1, double x2, double a, double b, double &y) override final
    {
        double xc1, xc2;
        p2_solve(a, b-a*beta*(alpha - 2), b*beta*(1.0-alpha),xc1, xc2);
        double y1 = pdf(xc1)/(a*xc1+b);
        double y2 = pdf(xc2)/(a*xc2+b);
        if((x1 < xc1 && xc1 < x2) && (x1 < xc2 && xc2 < x2)){
            y = std::min(y1, y2);
            return std::max<double>(y1, y2);
        } else if (x1 < xc1 && xc1 < x2){
            y = y1;
            return y1;
        } else if (x1 < xc2 && xc2 < x2){
            y = y2;
            return y2;
        } else{
            double xc = (x2+x1)*0.5;
            y = pdf(xc)/(a*xc+b);
            return y;
        }//*/
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

    virtual double left_aprox_pdf(double x) override final
    {
        return pdf(x);
    }
    virtual double right_aprox_pdf(double x) override final
    {
        return pdf(x);
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

    virtual double pdf(double x) override final
    {
        return nGammaDistributionTools::gamma_pdf_fast(x, alpha, beta);
    }

    virtual double qf(double x) override final
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
class GammaDistributionSplited final : public SplitedDistribution<GenType> 
{
private:
    double alpha;
    double beta;
    GammaDistribution<GenType> base_generator;

    GenType::GeneratorStateType init(double alpha, double beta, GenType::GeneratorStateType gs)
    {
        this->alpha = alpha;
        this->beta = beta;
        gs = this->init_splited(gs);
    }

protected:
    virtual double base_pdf(double x) override final
    {
        return nGammaDistributionTools::gamma_pdf_fast(x, alpha, beta);
    }
    virtual double base_cdf(double x) override final
    {
        return nGammaDistributionTools::gamma_cdf(x, alpha, beta);
    }

    virtual double base_qf(double q) override final
    {
        return nGammaDistributionTools::gamma_qf(q, alpha, beta);
    }

    virtual double base_gen() override final
    {
        return base_generator();
    }
    virtual GenType::GeneratorStateType init_base(GenType::GeneratorStateType gs) override final
    {
        base_generator =  GammaDistribution<GenType>(alpha, beta, gs);
        return gs;
    }

    virtual double computate_p(void)
    {
        if(alpha > 1.0){
            return nGammaDistributionTools::computate_p(alpha);
        } else{
            if(alpha <= 1.99427844665415544839844742611e-20){
                return 0;
            }
            //use exponential-logarithmic form due to fast grown of gamma(alpha)
            return exp(alpha*log(alpha)+1.0-alpha-lgamma(alpha));
        }
    }
    virtual double computate_s(void)
    {
        if(alpha <= 1.0){
            return 0;
        } else{
            return beta*(alpha-1.0);
        }
    }
    virtual double computate_lambda(void)
    {
        if(alpha < 1.0){
            return 1.0/(alpha*beta);
        }else{
            return 1.0/beta;
        }
    }

public:
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
};

//*/
}