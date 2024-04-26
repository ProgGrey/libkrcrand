#include "../../headers/generators.hpp"
#include "../../headers/gammaDistribution.hpp"
#include "../../headers/math.hpp"
#include "../3rdparty/cdflib.hpp"
#include <deque>
#include <algorithm>

using namespace krcrand;
using namespace std;


double nGammaDistributionTools::computate_p(double alpha)
{
    if(alpha <= 171.0){
        double pw = pow(alpha-1, alpha*0.5-0.5);
        return pw/tgamma(alpha)/exp(alpha-1)*pw;
    } else if(alpha > 137438953472.0){
       return 0.0;
    }else{
        long double a = static_cast<long double>(alpha);
        return static_cast<double>(expl((a-1.0)*logl(a-1.0)-lgammal(a)-a+1.0));
    }
}

double nGammaDistributionTools::gamma_cdf(double x, double alpha, double beta)
{
    int which = GAMMA_CDF;
    double p, q;
    double shape = alpha, rate = 1.0/beta;
    int status;
    double bound;
    n3rdPartyCDF::cdfgam(&which, &p, &q, &x, &shape, &rate, &status, &bound);
    if(status == 0){
        return p;
    } else{
        return -INFINITY;
    }
}

double nGammaDistributionTools::gamma_qf(double p, double alpha, double beta)
{
    int which = GAMMA_QF;
    double q = 1.0 - p;
    double shape = alpha, rate = 1.0/beta;
    int status;
    double bound;
    double x;
    n3rdPartyCDF::cdfgam(&which, &p, &q, &x, &shape, &rate, &status, &bound);
    if(status == 0){
        return x;
    } else{
        return -INFINITY;
    }
}

double nGammaDistributionTools::gamma_pdf(double _x_, double alpha, double beta)
{
    long double a = alpha;
    long double b = beta;
    long double x = _x_;
    if(x == 0.0){
        return (alpha < 1.0) ? INFINITY : ((alpha == 1.0) ? (1/beta) : 0);
    }
    return expl((a-1.0)*logl(x)-x/b-lgammal(a)-a*logl(b));
}

double nGammaDistributionTools::gamma_pdf_fast(double x, double alpha, double beta)
{
    return exp((alpha-1.0)*log(x)-x/beta-lgamma(alpha)-alpha*log(beta));
}