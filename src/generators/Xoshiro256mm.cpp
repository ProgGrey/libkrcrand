#include "../tools.hpp"
#include "../../headers/generators.hpp"
#include <cstring>

using namespace krcrand;

void Xoshiro256mmState::seed(uint64_t s)
{
    split_mix64_fill(s, state, 4);
}


void Xoshiro256mmState::next_state(uint64_t buf[4]) const
{
    uint64_t tmp = buf[1] << 17;

	buf[2] ^= buf[0];
	buf[3] ^= buf[1];
	buf[1] ^= buf[2];
	buf[0] ^= buf[3];

	buf[2] ^= tmp;

	buf[3] = rotl(buf[3], 45);
}

Xoshiro256mmState Xoshiro256mmState::jump(void) const
{
    static const uint64_t Jump[4] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
    uint64_t buf[State_Size];
    memcpy(buf, state, State_Size*sizeof(uint64_t));
    Xoshiro256mmState ret;
    ret.state[0] = 0;
    ret.state[1] = 0;
    ret.state[2] = 0;
    ret.state[3] = 0;
	for(unsigned int k = 0; k < sizeof(Jump) / sizeof(*Jump); k++){
		for(unsigned int j = 0; j < 64; j++) {
			if (Jump[k] & UINT64_C(1) << j) {
				ret.state[0] ^= buf[0];
				ret.state[1] ^= buf[1];
				ret.state[2] ^= buf[2];
				ret.state[3] ^= buf[3];
			}
			next_state(buf);	
		}
    }
    return ret;
}

const uint64_t* Xoshiro256mmState::raw(void) const
{
    return state;
}

Xoshiro256mmState Xoshiro256mmUniversal::set_state(const Xoshiro256mmState &state)
{
    memcpy(this->state, state.raw(), sizeof(uint64_t)*state.State_Size);
    return state.jump();
}

void Xoshiro256mmUniversal::fill(void)
{
    for(uint_fast16_t k = 0; k < Generator_Buff_Size; k++){
        this->buffer[k] = gen();
    }
}

uint64_t Xoshiro256mmUniversal::gen(void)
{
    uint64_t value = rotl(state[1] * 5, 7) * 9;

	uint64_t tmp = state[1] << 17;

	state[2] ^= state[0];
	state[3] ^= state[1];
	state[1] ^= state[2];
	state[0] ^= state[3];

	state[2] ^= tmp;

	state[3] = rotl(state[3], 45);

	return value;
}

#ifdef LIBKRCRAND_ENABLE_SSE2

Xoshiro256mmState Xoshiro256mmSSE2::set_state(const Xoshiro256mmState &state)
{
	uint64_t s1[state.State_Size];
	uint64_t s2[state.State_Size];
	uint64_t buf[2];
    memcpy(s1, state.raw(), sizeof(uint64_t)*state.State_Size);
	auto tmp = state.jump();
	memcpy(s2, tmp.raw(), sizeof(uint64_t)*tmp.State_Size);
	for(unsigned int k = 0; k < state.State_Size; k++){
		buf[0] = s1[k];
		buf[1] = s2[k];
		memcpy(&(this->state[k]), buf, sizeof(uint64_t)*2);
	}
    return tmp.jump();
}

__m128i Xoshiro256mmSSE2::gen(void)
{
	//__m128i value = rotl_sse2(mull_u64_sse2(state[1], _mm_set1_epi64x(5)), 7);
	//value = mull_u64_sse2(value, _mm_set1_epi64x(9));
	__m128i value = _mm_add_epi64(state[1], _mm_slli_epi64(state[1], 2));//Equals multiplying by 5
	value = rotl_sse2(value, 7);
	value = _mm_add_epi64(value, _mm_slli_epi64(value, 3));//Equals multiplying by 9
	__m128i tmp = _mm_slli_epi64(state[1], 17);
	
	state[2] = _mm_xor_si128(state[2], state[0]);
	state[3] = _mm_xor_si128(state[3], state[1]);
	state[1] = _mm_xor_si128(state[1], state[2]);
	state[0] = _mm_xor_si128(state[0], state[3]);

	state[2] = _mm_xor_si128(state[2], tmp);

	state[3] = rotl_sse2(state[3], 45);
	return value;
}

void Xoshiro256mmSSE2::fill(void)
{
	for(unsigned int k = 0; k < Generator_Buff_Size; k+=2){
		__m128i tmp = gen();
		memcpy(buffer + k, &tmp, sizeof(uint64_t)*2);
	}
	if(Generator_Buff_Size & 1){
		__m128i tmp = gen();
		memcpy(buffer + Generator_Buff_Size - 1, &tmp, sizeof(uint64_t));
	}
	//*/
}

#endif



#ifdef LIBKRCRAND_ENABLE_AVX2

Xoshiro256mmState Xoshiro256mmAVX2::set_state(const Xoshiro256mmState &state)
{
	uint64_t s[4][state.State_Size];
	uint64_t buf[4];
	Xoshiro256mmState tmp = state;
	for(unsigned int k =0; k < 4; k++){
		memcpy(s[k], tmp.raw(), sizeof(uint64_t)*tmp.State_Size);
		tmp = tmp.jump();
	}
	for(unsigned int k = 0; k < state.State_Size; k++){
		for(unsigned int i = 0; i < 4; i++){
			buf[i] = s[i][k];
		}
		memcpy(&(this->state[k]), buf, sizeof(uint64_t)*4);
	}
    return tmp;
}

__m256i Xoshiro256mmAVX2::gen(void)
{
	__m256i value = _mm256_add_epi64(state[1], _mm256_slli_epi64(state[1], 2));//Equals multiplying by 5
	value = rotl_avx2(value, 7);
	value = _mm256_add_epi64(value, _mm256_slli_epi64(value, 3));//Equals multiplying by 9
	__m256i tmp = _mm256_slli_epi64(state[1], 17);
	
	state[2] = _mm256_xor_si256(state[2], state[0]);
	state[3] = _mm256_xor_si256(state[3], state[1]);
	state[1] = _mm256_xor_si256(state[1], state[2]);
	state[0] = _mm256_xor_si256(state[0], state[3]);

	state[2] = _mm256_xor_si256(state[2], tmp);

	state[3] = rotl_avx2(state[3], 45);
	return value;
}

void Xoshiro256mmAVX2::fill(void)
{
	for(unsigned int k = 0; k < Generator_Buff_Size; k+=4){
		__m256i tmp = gen();
		memcpy(buffer + k, &tmp, sizeof(uint64_t)*4);
	}
	if((Generator_Buff_Size % 4) != 0){
		__m256i tmp = gen();
		memcpy(buffer + Generator_Buff_Size - (Generator_Buff_Size % 4) , &tmp, (Generator_Buff_Size % 4)*sizeof(uint64_t));
	}
	//*/
}

#endif

#ifdef LIBKRCRAND_ENABLE_AVX512F

Xoshiro256mmState Xoshiro256mmAVX512F::set_state(const Xoshiro256mmState &state)
{
	uint64_t s[8][state.State_Size];
	uint64_t buf[8];
	Xoshiro256mmState tmp = state;
	for(unsigned int k =0; k < 8; k++){
		memcpy(s[k], tmp.raw(), sizeof(uint64_t)*tmp.State_Size);
		tmp = tmp.jump();
	}
	for(unsigned int k = 0; k < state.State_Size; k++){
		for(unsigned int i = 0; i < 8; i++){
			buf[i] = s[i][k];
		}
		memcpy(&(this->state[k]), buf, sizeof(uint64_t)*8);
	}
    return tmp;
}

__m512i Xoshiro256mmAVX512F::gen(void)
{
	__m512i value = _mm512_add_epi64(state[1], _mm512_slli_epi64(state[1], 2));//Equals multiplying by 5
	value = rotl_avx512f(value, 7);
	value = _mm512_add_epi64(value, _mm512_slli_epi64(value, 3));//Equals multiplying by 9
	__m512i tmp = _mm512_slli_epi64(state[1], 17);
	
	state[2] = _mm512_xor_si512(state[2], state[0]);
	state[3] = _mm512_xor_si512(state[3], state[1]);
	state[1] = _mm512_xor_si512(state[1], state[2]);
	state[0] = _mm512_xor_si512(state[0], state[3]);

	state[2] = _mm512_xor_si512(state[2], tmp);

	state[3] = rotl_avx512f(state[3], 45);
	return value;
}

void Xoshiro256mmAVX512F::fill(void)
{
	for(unsigned int k = 0; k < Generator_Buff_Size; k+=8){
		__m512i tmp = gen();
		memcpy(buffer + k, &tmp, sizeof(uint64_t)*8);
		//_mm512_store_epi64(buffer + k, tmp);
	}
	if((Generator_Buff_Size % 8) != 0){
		__m512i tmp = gen();
		memcpy(buffer + Generator_Buff_Size - (Generator_Buff_Size % 8) , &tmp, (Generator_Buff_Size % 8)*sizeof(uint64_t));
	}
	//*/
}

#endif