#include "../headers/generators.hpp"
#include <iostream>
#include <random>

using namespace std;
using namespace krcrand;

int main()
{
    random_device rd;
    Xoshiro256mmState st(rd());
    Xoshiro256mm gen(st);
    while(true){
        uint64_t val = gen();
        cout.write(reinterpret_cast<const char*>(&val), sizeof(uint64_t));   
    }
    return 0;
}