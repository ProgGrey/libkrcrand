#include "platform.hpp"
#include "math.hpp"
#include <typeinfo>

namespace krcrand{

template<typename GenType, int is_inverted> class ExponentialDistribution
{
private:
    unsigned int pos;
    DECL_KRCRAND_ALIGN double buf[Generator_Buff_Size];
    double mult;
    auto entropy()
    {
        if(is_inverted){
            return uniform01_exclude0(generator.gen_n());
        } else {
            return uniform01_exclude0(generator.gen());
        }
    }
    auto gen()
    {
        return type_depend_mull(unsafe_log(entropy()), mult);
    }
    void fill(void)
    {
        unsigned int step = sizeof(std::invoke_result_t<decltype(&ExponentialDistribution<GenType, is_inverted>::gen), ExponentialDistribution<GenType, is_inverted>>) / sizeof(double);
        for(unsigned int k = 0; k < Generator_Buff_Size; k+= step){
            auto val = gen();
            memcpy(buf + k, &val, sizeof(val));
        }
    }
public:
    GenType generator;
    double operator()(){
        pos++;
        if(pos>=Generator_Buff_Size){
            pos = 0;
            fill();
        }
        return(buf[pos]);
    }
    explicit ExponentialDistribution(double lambda){
        pos = Generator_Buff_Size;
        generator.seed(0);
        mult = -1.0/lambda;
    }
    ExponentialDistribution(uint64_t seed, double lambda){
        pos = Generator_Buff_Size;
        generator.seed(seed);
        mult = -1.0/lambda;
    }
};

}