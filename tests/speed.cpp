#include "../headers/libkrcrand.hpp"
#include <iostream>
#include <chrono>
#include <random>

using namespace std;
using namespace krcrand;

int main()
{
    random_device rd;
    mt19937_64 gen0(rd());
    Xoshiro256mmState st(rd());
    Xoshiro256mmUniversal gen1(st);
    Xoshiro256mmSSE2 gen2(st);
    double res = 0;
    uint64_t N = 10000000000;
    auto start = chrono::system_clock::now();
    for(uint64_t k = 0; k < N; k++){
        res += (double)gen0()/N;
    }
    auto end = chrono::system_clock::now();
    auto diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "mt19937_64:" << diff.count() << "ms\t" << res << endl;
    start = chrono::system_clock::now();
    for(uint64_t k = 0; k < N; k++){
        res += (double)gen1()/N;
    }
    end = chrono::system_clock::now();
    diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Xoshiro256mmUniversal:" << diff.count() << "ms\t" << res << endl;
    start = chrono::system_clock::now();
    for(uint64_t k = 0; k < N; k++){
        res += (double)gen2()/N;
    }
    end = chrono::system_clock::now();
    diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Xoshiro256mmSSE2:" << diff.count() << "ms\t" << res << endl;
    return 0;
}