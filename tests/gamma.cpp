#include "../headers/libkrcrand.hpp"
#include "../headers/generators.hpp"
#include "../headers/gammaDistribution.hpp"
#include <iostream>
#include <format>

using namespace std;
using namespace krcrand;

int main()
{
    GammaDistributionSplited<Xoshiro256mm> gen(3,2);
    for(unsigned int j =0; j < 1000; j++){
        for(unsigned int k =0; k < 100; k++){
            cout << format("{}", gen()) << ", ";
        }
        cout << endl;
    }
    cout << endl;
    cout << gen.max_fails  << endl;
    cout << gen.sum_fails << endl;
    return 0;   
}