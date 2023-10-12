#include "../tools.hpp"
#include "../../headers/generators.hpp"
#include <cstring>

void Xoshiro256mmState::seed(uint64_t s)
{
    split_mix64_fill(s, state, 4);
}


void Xoshiro256mmState::next_state(uint64_t buf[4])
{
    uint64_t tmp = buf[1] << 17;

	buf[2] ^= buf[0];
	buf[3] ^= buf[1];
	buf[1] ^= buf[2];
	buf[0] ^= buf[3];

	buf[2] ^= tmp;

	buf[3] = rotl(buf[3], 45);
}

Xoshiro256mmState Xoshiro256mmState::jump(void)
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

const uint64_t* Xoshiro256mmState::raw(void)
{
    return state;
}

Xoshiro256mmState Xoshiro256mmUniversal::set_state(Xoshiro256mmState &state)
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