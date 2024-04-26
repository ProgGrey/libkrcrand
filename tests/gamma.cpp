//#include "../headers/libkrcrand.hpp"
#include <iostream>
#include "../headers/generators.hpp"
#include "../headers/gammaDistribution.hpp"
#include <format>
#include <vector>
#include <deque>
#include <algorithm>
#include <iomanip>
#include <quadmath.h>

using namespace std;
using namespace krcrand;

double horner1_(double x, const double* coefs, uint_fast8_t n)
{
    double val = x + coefs[0];
    for (uint_fast8_t k = 1; k < n; k++) {
        val = fma(x, val, coefs[k]);
    }
    return val;
}

double lambertw1(double x)
{
    //static const double[10][7] ={{}};
    return x + 1;
}


double computate_p(double alpha)
{
    if(alpha <= 171.0){
        double pw = pow(alpha-1, alpha*0.5-0.5);
        return pw/tgamma(alpha)/exp(alpha-1)*pw;
    } else if(alpha > 137438953472.0){
       return 0.0;
    }else{
        double a = alpha;
        return exp((a-1.0)*log(a-1.0)-lgamma(a)-a+1.0);
    }
}
// 
int main()
{
    std::cout << std::setprecision (std::numeric_limits<double>::digits10 + 1);
    Xoshiro256mm::GeneratorStateType inter(1);
    //GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, inter);
    GammaDistribution<Xoshiro256mm> gen(0.1, 0.2, inter);
    for(unsigned int j =0; j < 100; j++){
        for(unsigned int k =0; k < 10; k++){
            cout << gen() << ", ";
            //gen();
        }
        cout << endl;
    }
    /*
    GammaDistributionSplited<Xoshiro256mm> gen2(nextafter(1,2), 1, inter);
    
    cout << gen.max_fails << endl;
    cout << endl;
    for(unsigned int j =0; j < 1; j++){
        for(unsigned int k =0; k < 10; k++){
            cout << format("{}", gen2()) << ", ";
        }
        cout << endl;
    }
    cout << endl;
    cout << gen.max_fails  << endl;
    cout << gen.sum_fails << endl;
    */
    return 0;   
}