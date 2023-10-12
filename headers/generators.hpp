#ifndef __LIBKRCRAND_GENERATORS_HPP__
#define __LIBKRCRAND_GENERATORS_HPP__

#include <cstdint>
#include "constants.hpp"
#include "platform.hpp"

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#endif

namespace krcrand{

class InternalState
{
private:
    
public:
    virtual void seed(uint64_t s) = 0;
    virtual const uint64_t* raw(void) = 0;
};

template<typename InternalStateClass> class Generator
{
protected:
    // Buffer with random numbers
    uint64_t buffer[Generator_Buff_Size];
    uint16_t position;
    // Fill buffer with random numbers
    virtual void fill(void) = 0;
public:
    // Return next random value
    inline uint64_t operator()(void)
    {
        if(position >= Generator_Buff_Size){
            fill();
            position = 0;
        }
        uint64_t value = buffer[position];
        position++;
        return value;
    }

    // Returns state for next generator for chained initialization.
    virtual InternalStateClass set_state(InternalStateClass &state) = 0;

    Generator()
    {
        position = Generator_Buff_Size;
    }
};

// Xoshiro256** interal state class
class Xoshiro256mmState : public InternalState
{
public:
    // Size of internal state space
    static const unsigned int State_Size = 4;
    // Equals 2^128 number generations
    Xoshiro256mmState jump(void);
    // Set seed
    virtual void seed(uint64_t s) override final;
    // Return pointer to internal state
    virtual const uint64_t* raw(void) override final;
private:
    void next_state(uint64_t buf[State_Size]);
    uint64_t state[State_Size];
};

// Xoshiro256** universal version for all platforms
class Xoshiro256mmUniversal : public Generator<Xoshiro256mmState>
{
private:
    uint64_t state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
    uint64_t gen(void);
public:
    virtual Xoshiro256mmState set_state(Xoshiro256mmState &state) override final;
};

#ifdef LIBKRCRAND_ENABLE_SSE2
// Xoshiro256** SSE2 version
class Xoshiro256mmSSE2 : public Generator<Xoshiro256mmState>
{
private:
    __m128i  state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
    __m128i gen(void);
public:
    virtual Xoshiro256mmState set_state(Xoshiro256mmState &state) override final;
};
#endif

}
#endif