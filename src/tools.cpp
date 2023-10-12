#include "tools.hpp"

void split_mix64_fill(uint64_t seed, uint64_t *buf, unsigned int size)
{
    for(unsigned int k = 0; k < size; k++){
        uint64_t val = (seed += 0x9e3779b97f4a7c15);
        val = (val ^ (val >> 30)) * 0xbf58476d1ce4e5b9;
        val = (val ^ (val >> 27)) * 0x94d049bb133111eb;
        buf[k] =  val ^ (val >> 31);
    }
}

uint64_t rotl(const uint64_t x, uint_fast8_t k) 
{
	return (x << k) | (x >> (64 - k));
}
