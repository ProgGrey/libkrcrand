#include "../headers/libkrcrand.hpp"
#include "../headers/generators.hpp"
#include "../headers/gammaDistribution.hpp"
#include <iostream>
#include <format>
#include <vector>
#include <deque>
#include <algorithm>

using namespace std;
using namespace krcrand;


int main()
{
    Xoshiro256mm::GeneratorStateType inter(1);
    GammaDistributionSplited<Xoshiro256mm> gen(nextafter(1,2), 1, inter);
    GammaDistributionSplited<Xoshiro256mm> gen2(nextafter(1,2), 1, inter);
    for(unsigned int j =0; j < 1; j++){
        for(unsigned int k =0; k < 10; k++){
            cout << format("{}", gen()) << ", ";
        }
        cout << endl;
    }
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
    return 0;   
}