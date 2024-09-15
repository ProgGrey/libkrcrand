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
    /*
    double x1,x2,x3;
    cout << sizeof(long double) << endl;
    std::cout << std::setprecision (std::numeric_limits<double>::digits10 + 1);
    cout << p3_solve_cardano(1e-15,-2,-10,0.5, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    cout << p3_solve_modern(4,-2,-10,0.5, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    cout << p3_solve(-1,1,1,0.5, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    cout << p3_solve(1e-15,-2,-10,0.5, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    cout << p3_solve(-1,2,1,-1, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    cout << p3_solve(-1,3.52,-2.7968,-0.541696, x1, x2, x3) << '\t' << x1 << '\t' << x2 << '\t' << x3 << endl;
    //*/
    Xoshiro256mm::GeneratorStateType inter(1);
    //GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, inter);
    GammaDistribution<Xoshiro256mm> gen(1.2, 10, inter);
    cout.precision(9);
    for(unsigned int j =0; j < 100; j++){
        for(unsigned int k =0; k < 10; k++){
            cout << gen() << ", ";
            //gen();
        }
        cout << endl;
    }
    cout << gen.counter << endl;
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