#include "platform.hpp"
#include <type_traits>
#include <typeinfo>
#include <iostream>

namespace krcrand{

template<typename GenType, int is_inverted> class ExponentialDistribution
{
private:
    GenType generator;
    unsigned int pos = Generator_Buff_Size;
    double buf[Generator_Buff_Size];
    double mult;
    auto entropy()
    {
        if(is_inverted){
            return generator.gen_n();
        } else {
            return generator.gen();
        }
    }
    auto gen()
    {
        return type_depend_mull(type_depend_log(entropy()), mult);
    }
    void fill(void)
    {
        if(typeid(typename std::invoke_result<GenType>::type) == typeid(uint64_t)){
            for(unsigned int k = 0; k < Generator_Buff_Size; k++){
                buf[k] = gen();
            }
        } else {
            throw;
        }
    }
public:
    double operator()(){
        pos++;
        if(pos>=Generator_Buff_Size){
            pos = 0;
            fill();
        }
        return(buf[pos]);
    }
    ExponentialDistribution(double lambda){
        mult = -1.0/lambda;
    }
};

}