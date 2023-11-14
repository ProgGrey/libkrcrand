#include "../headers/libkrcrand.hpp"
#include "../headers/exponentialDistribution.hpp"
#include <iostream>
#include <chrono>
#include <random>
#include <functional>

using namespace std;
using namespace krcrand;

#define test_generator(name, gen) \
    start = chrono::system_clock::now(); \
    for(uint64_t k = 0; k < N; k++){ \
        gen();\
    }\
    end = chrono::system_clock::now();\
    diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);\
    cout << name << diff.count() << "ms\t" << gen() << endl;

#define test_exp_dist(name, gen)\
    start = chrono::system_clock::now(); \
    mean = 0;\
    for(uint64_t k = 0; k < N; k++){ \
        mean += gen()/N;\
    }\
    end = chrono::system_clock::now();\
    diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);\
    cout << name << diff.count() << "ms\t" << 1/mean << endl;

int main()
{
    random_device rd;
    mt19937_64 gen0(rd());
    exponential_distribution<double> exp_dist(1.2);
    auto exp0 = bind(exp_dist, gen0);
    //Xoshiro256mmState st(rd());
    Xoshiro256mmState st(1);
    Xoshiro256mmUniversalStable gen1(st);
    
    double mean;
    uint64_t N = 1000000000;
    chrono::_V2::system_clock::time_point start;
    chrono::_V2::system_clock::time_point end;
    chrono::milliseconds diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    test_generator("mt19937_64: ", gen0)
    test_exp_dist("mt19937_64 Exp: ", exp0)
    test_generator("Xoshiro256mmUniversalStable: ", gen1)
    ExponentialDistribution<Xoshiro256mmUniversalStable, 0> exp1(1.2);
    exp1.generator.set_state(st);
    test_exp_dist("Xoshiro256mmUniversalStable Exp: ", exp1)
#ifdef LIBKRCRAND_ENABLE_SSE2
    Xoshiro256mmSSE2stable gen2(st);
    test_generator("Xoshiro256mmSSE2stable: ", gen2)
    ExponentialDistribution<Xoshiro256mmSSE2stable, 0> exp2(1.2);
    exp2.generator.set_state(st);
    test_exp_dist("Xoshiro256mmSSE2stable Exp: ", exp2)
#endif
#ifdef LIBKRCRAND_ENABLE_AVX2
    Xoshiro256mmAVX2stable gen3(st);
    test_generator("Xoshiro256mmAVX2stable: ", gen3)
    ExponentialDistribution<Xoshiro256mmAVX2stable, 0> exp3(1.2);
    exp3.generator.set_state(st);
    test_exp_dist("Xoshiro256mmAVX2stable Exp: ", exp3)
#endif
#ifdef LIBKRCRAND_ENABLE_AVX512F
    Xoshiro256mmAVX512Fstable gen4(st);
    test_generator("Xoshiro256mmAVX512Fstable: ", gen4)
    ExponentialDistribution<Xoshiro256mmAVX512Fstable, 0> exp4(1.2);
    exp4.generator.set_state(st);
    test_exp_dist("Xoshiro256mmAVX512Fstable Exp: ", exp4)
#endif
    Xoshiro256mm gen5(st);
    test_generator("Xoshiro256mm: ", gen5);
    ExponentialDistribution<Xoshiro256mm, 0> exp5(1.2);
    exp5.generator.set_state(st);
    test_exp_dist("Xoshiro256mm Exp: ", exp5)
    
    return 0;
}