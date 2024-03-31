#include "../../headers/generators.hpp"
#include "../../headers/gammaDistribution.hpp"
#include "../../headers/math.hpp"
#include <deque>
#include <algorithm>

using namespace krcrand;
using namespace std;

inline double stage1_fast(double alpha, uint_fast32_t ai)
{
    double res = 1.0;
    for(uint_fast32_t k = 1; k < ai/2; k++){
        res *= (alpha-1.0)*M_E_m2*(alpha-1.0)/((alpha-k)*(alpha-(ai-k)));
    }
    return res;
}

inline double stage1_slow(double alpha, uint_fast32_t ai)
{
    deque<double> deq1;
    for(uint_fast32_t k = 1; k < ai/2; k++){
        deq1.push_back((alpha-1.0)*M_E_m2*(alpha-1.0)/((alpha-k)*(alpha-(ai-k))));
    }
    sort(deq1.begin(), deq1.end());
    while(deq1.size() > 1){
        deque<double> deq2;
        while(deq1.size() >= 2){
            deq2.push_back(deq1.front()*deq1.back());
            deq1.pop_front();
            deq1.pop_back();
        }
        if(deq1.size() == 1){
            deq2.push_back(deq1.front());
        }
        deq1 = deq2;
        sort(deq1.begin(), deq1.end());

    };
    return deq1.front();
}

double nGammaDistributionTools::computate_p(double alpha)
{
    if(alpha <= 171.0){
        double pw = pow(alpha-1, alpha*0.5-0.5);
        return pw/tgamma(alpha)/exp(alpha-1)*pw;
    } else if(alpha > 4294967295.0){
        return 0;
    }
    uint_fast32_t  ai = static_cast<uint_fast32_t>(alpha);

    double res;
    if(alpha < 4867){
        res = stage1_fast(alpha, ai);
    } else{
        res = stage1_slow(alpha, ai);
    }
    
    if(ai & 1){
        res *= (alpha-1.0)*M_E_m2/((alpha-ai/2-1.0))*(alpha-1.0)/((alpha-ai/2));
    } else{
        res *= (alpha-1.0)/((alpha-ai/2))*M_E_m1;
    }
    return pow(alpha-1.0, alpha-ai)*res/(exp(alpha-ai)*tgamma(alpha-ai+1.0));
}