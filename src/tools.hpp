#ifndef __LIBKRCRAND_TOOLS_HPP__
#define __LIBKRCRAND_TOOLS_HPP__

#include <cstdint>


void split_mix64_fill(uint64_t seed, uint64_t *buf, unsigned int  size);
uint64_t rotl(const uint64_t x, uint_fast8_t k);

#endif

