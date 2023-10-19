#include "../headers/libkrcrand.hpp"
#include <iostream>
#include <chrono>
#include <random>

using namespace std;
using namespace krcrand;

#define test_generator(name, gen) \
    start = chrono::system_clock::now(); \
    for(uint64_t k = 0; k < N; k++){ \
        gen();\
    }\
    end = chrono::system_clock::now();\
    diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);\
    cout << name << diff.count() << "ms\t" << gen() << endl;\

int main()
{
    random_device rd;
    mt19937_64 gen0(rd());
    Xoshiro256mmState st(rd());
    Xoshiro256mmUniversalStable gen1(st);

    uint64_t N = 1000000000;
    auto start = chrono::system_clock::now();
    auto end = chrono::system_clock::now();
    auto diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    test_generator("mt19937_64: ", gen0)
    test_generator("Xoshiro256mmUniversalStable: ", gen1)
#ifdef LIBKRCRAND_ENABLE_SSE2
    Xoshiro256mmSSE2stable gen2(st);
    test_generator("Xoshiro256mmSSE2stable: ", gen2)
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    Xoshiro256mmAVX2stable gen3(st);
    test_generator("Xoshiro256mmAVX2stable: ", gen3)
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    Xoshiro256mmAVX512F gen4(st);
#endif
    return 0;
}