#ifndef __LIBKRCRAND_GENERATORS_HPP__
#define __LIBKRCRAND_GENERATORS_HPP__

#include <cstdint>
#include <cstring>
#include "constants.hpp"
#include "platform.hpp"

#ifdef LIBKRCRAND_ENABLE_SSE2
#include <emmintrin.h>
#endif

#if LIBKRCRAND_ENABLE_AVX2 || LIBKRCRAND_ENABLE_AVX512F
#include <immintrin.h>
#endif

#include <iostream>

namespace krcrand{

class InternalState
{
private:
    
public:
    virtual void seed(uint64_t s) = 0;
    virtual const uint64_t* raw(void) const = 0;
};

template<typename InternalStateClass> class Generator
{
protected:
    // Buffer with random numbers
    DECL_KRCRAND_ALIGN uint64_t buffer[Generator_Buff_Size];
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
    virtual InternalStateClass set_state(const InternalStateClass &state) = 0;

    Generator()
    {
        position = Generator_Buff_Size;
        memset(buffer, 0, Generator_Buff_Size*sizeof(uint64_t));
    }
};

// Xoshiro256** interal state class
class Xoshiro256mmState : public InternalState
{
public:
    // Size of internal state space
    static const unsigned int State_Size = 4;
    // Equals 2^128 number generations
    Xoshiro256mmState jump(void) const;
    // Set seed
    virtual void seed(uint64_t s) override final;
    // Return pointer to internal state
    virtual const uint64_t* raw(void) const override final;
    explicit Xoshiro256mmState(uint64_t val)
    {
        seed(val);
    }
    Xoshiro256mmState(){}
private:
    void next_state(uint64_t buf[State_Size]) const;
    uint64_t state[State_Size];
};

// Xoshiro256** universal version for all platforms
class Xoshiro256mmUniversal : public Generator<Xoshiro256mmState>
{
private:
    uint64_t state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
public:
    uint64_t gen(void);
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final;
    Xoshiro256mmUniversal(){}
    explicit Xoshiro256mmUniversal(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};

class Xoshiro256mmUniversalStable: public Generator<Xoshiro256mmState>
{
    private:
    uint8_t pos = 0;
    Xoshiro256mmUniversal gens[COMPAB_COUNT_UNIV];
    virtual void fill(void) override final
    {
        uint_fast8_t gnum = 0;
        for(unsigned int k = 0; k < Generator_Buff_Size; k++){
            buffer[k] = gens[gnum]();
            //std::cout << (int)gnum << '\t';
            gnum++;
            if(gnum >= COMPAB_COUNT_UNIV){
                gnum = 0;
            }
            //buffer[k] = 0;
        }
    }
    public:
    uint64_t gen(void)
    {
        uint64_t ret = gens[pos].gen();
        pos++;
        if(pos >= COMPAB_COUNT_UNIV){
            pos = 0;
        }
        return ret;
    }
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final
    {
        Xoshiro256mmState st = state;
        for(unsigned int k = 0; k < COMPAB_COUNT_UNIV; k++){
            gens[k].set_state(st);
            st = st.jump();
        }
        return st;
    }
    Xoshiro256mmUniversalStable(){}
    explicit Xoshiro256mmUniversalStable(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};

#ifdef LIBKRCRAND_ENABLE_SSE2
// Xoshiro256** SSE2 version
class Xoshiro256mmSSE2 : public Generator<Xoshiro256mmState>
{
private:
    alignas(16) __m128i  state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
public:
    __m128i gen(void);
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final;
    Xoshiro256mmSSE2(){}
    explicit Xoshiro256mmSSE2(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};

class Xoshiro256mmSSE2stable: public Generator<Xoshiro256mmState>
{
    private:
    uint8_t pos = 0;
    Xoshiro256mmSSE2 gens[COMPAB_COUNT_SSE2];
    virtual void fill(void) override final
    {
        uint_fast8_t gnum = 0;
        for(unsigned int k = 0; k < Generator_Buff_Size; k+=2){
            for(unsigned int j = 0; j < 2; j++){
                buffer[k+j] = gens[gnum]();
            }
            gnum++;
            if(gnum >= COMPAB_COUNT_SSE2){
                gnum = 0;
            }
        }
    }
    public:
    __m128i gen(void)
    {
        __m128i ret = gens[pos].gen();
        pos++;
        if(pos >= COMPAB_COUNT_AVX2){
            pos = 0;
        }
        return ret;
    }
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final
    {
        Xoshiro256mmState st = state;
        for(unsigned int k = 0; k < COMPAB_COUNT_SSE2; k++){
            st = gens[k].set_state(st);
        }
        return st;
    }
    Xoshiro256mmSSE2stable(){}
    explicit Xoshiro256mmSSE2stable(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};
#endif

#ifdef LIBKRCRAND_ENABLE_AVX2
// Xoshiro256** AVX2 version
class Xoshiro256mmAVX2 : public Generator<Xoshiro256mmState>
{
private:
    alignas(32) __m256i  state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
public:
    __m256i gen(void);
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final;
    Xoshiro256mmAVX2(){}
    explicit Xoshiro256mmAVX2(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};

class Xoshiro256mmAVX2stable: public Generator<Xoshiro256mmState>
{
    private:
    uint8_t pos = 0;
    Xoshiro256mmAVX2 gens[COMPAB_COUNT_AVX2];
    virtual void fill(void) override final
    {
        uint_fast8_t gnum = 0;
        for(unsigned int k = 0; k < Generator_Buff_Size; k+=4){
            for(unsigned int j = 0; j < 4; j++){
                buffer[k+j] = gens[gnum]();
            }
            gnum++;
            if(gnum >= COMPAB_COUNT_AVX2){
                gnum = 0;
            }
        }
    }
    public:
    __m256i gen(void)
    {
        __m256i ret = gens[pos].gen();
        pos++;
        if(pos >= COMPAB_COUNT_AVX2){
            pos = 0;
        }
        return ret;
    }
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final
    {
        Xoshiro256mmState st = state;
        for(unsigned int k = 0; k < COMPAB_COUNT_AVX2; k++){
            st = gens[k].set_state(st);
        }
        return st;
    }
    Xoshiro256mmAVX2stable(){}
    explicit Xoshiro256mmAVX2stable(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
// Xoshiro256** AVX2 version
class Xoshiro256mmAVX512F: public Generator<Xoshiro256mmState>
{
private:
    alignas(64) __m512i  state[Xoshiro256mmState::State_Size];
    virtual void fill(void) override final;
public:
    __m512i gen(void);
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final;
    Xoshiro256mmAVX512F(){}
    explicit Xoshiro256mmAVX512F(const Xoshiro256mmState &state) 
    {
        set_state(state);
    }
};

#if COMPAB_COUNT_AVX512F == 1
typedef Xoshiro256mmAVX512F Xoshiro256mmAVX512Fstable;
#else
class Xoshiro256mmAVX512Fstable: public Generator<Xoshiro256mmState>
{
    private:
    uint8_t pos = 0;
    Xoshiro256mmAVX512F gens[COMPAB_COUNT_AVX512F];
    virtual void fill(void) override final
    {
        uint_fast8_t gnum = 0;
        for(unsigned int k = 0; k < Generator_Buff_Size; k+=8){
            for(unsigned int j = 0; j < 8; j++){
                buffer[k+j] = gens[gnum]();
            }
            gnum++;
            if(gnum >= COMPAB_COUNT_AVX512F){
                gnum = 0;
            }
        }
    }
    public:
    __m512i gen(void)
    {
        __m512i ret = gens[pos].gen();
        pos++;
        if(pos >= COMPAB_COUNT_AVX512F){
            pos = 0;
        }
        return ret;
    }
    virtual Xoshiro256mmState set_state(const Xoshiro256mmState &state) override final
    {
        Xoshiro256mmState st = state;
        for(unsigned int k = 0; k < COMPAB_COUNT_AVX512F; k++){
            st = gens[k].set_state(st);
        }
        return st;
    }
    Xoshiro256mmAVX512Fstable(){}
    explicit Xoshiro256mmAVX512Fstable(const Xoshiro256mmState &state)
    {
        set_state(state);
    }
};
#endif
#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F
typedef  Xoshiro256mmAVX512Fstable Xoshiro256mm;
#elif defined(LIBKRCRAND_ENABLE_AVX2)
typedef  Xoshiro256mmAVX2stable Xoshiro256mm;
#elif defined(LIBKRCRAND_ENABLE_SSE2)
typedef  Xoshiro256mmSSE2stable Xoshiro256mm;
#else
typedef  Xoshiro256mmUniversalStable Xoshiro256mm;
#endif

}
#endif