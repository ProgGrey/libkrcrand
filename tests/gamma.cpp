#include "../headers/libkrcrand.hpp"
#include "../headers/generators.hpp"
#include "../headers/gammaDistribution.hpp"
#include <iostream>
#include <format>

using namespace std;
using namespace krcrand;

int main()
{
    GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, 1);
    for(unsigned int j =0; j < 1; j++){
        for(unsigned int k =0; k < 10; k++){
            cout << format("{}", gen()) << ", ";
        }
        cout << endl;
    }
    cout << endl;
    cout << gen.max_fails  << endl;
    cout << gen.sum_fails << endl;
    return 0;   
}