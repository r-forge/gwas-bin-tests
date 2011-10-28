/*
 * random.hpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#ifndef RANDOM_HPP_
#define RANDOM_HPP_


#include <time.h>
#include <sfmt.h>


// use rand() : not multithread-safe, uses global status
class Random {
public:
	Random(int seed=time(NULL), int seed2 = 0) { srand(seed); }
	~Random() {}

	//inline double random() { return rand()/(RAND_MAX+1.0); }


public: // ========== CLASS METHODS ==========
	static inline double random() { return rand()/(RAND_MAX+1.0); }
	static inline int random_int_below(int n) {
		return (int) ( ((double)n)* ( rand()/(RAND_MAX+1.0) ) );
	}
};

//  use random_r() : nmultithread-safe
class Random_R {
public:
	Random_R(int seed=time(NULL), int seed2 = 0) {
		_random_data.state = 0;
		initstate_r(seed, _statebuf, state_length, &_random_data);
	}
	~Random_R() {}

	inline double random() {
		int32_t res;
		random_r(&_random_data, &res);
		return res / (RAND_MAX+1.0);
	}
	inline int random_int_below(int n) {
		int32_t res;
		random_r(&_random_data, &res);
		return (int) ( ((double)n)* ( res/(RAND_MAX+1.0) ) );
	}
private:
	static const int state_length = 256;
	struct random_data _random_data;
	char _statebuf[state_length];
};

//class CRandomSFMT0_opt : public CRandomSFMT0 {
//public:
//	CRandomSFMT0_opt(unsigned int seed) : CRandomSFMT0(seed) {}
//	uint32_t BRandom();
//	int IRandom(int min, int max);
//	//uint32_t MotherBits();
//
////	void Generate();
//};
//
//class RandomSFMT{
//public:
//	RandomSFMT(unsigned int seed = time(NULL)) : _impl(seed) {}
//	inline double random() { return _impl.Random(); }
//	inline int random_int_below(int n) {
//		return _impl.IRandom(0, n-1);
//	}
//
//private:
//	CRandomSFMT0_opt _impl;
//};


class CRandomSFMT0_karl : public CRandomSFMT0 {
public:
	CRandomSFMT0_karl(int seed) : CRandomSFMT0(seed) {}
	CRandomSFMT0_karl(int seed1, int seed2);
	uint32_t BRandom();

	unsigned int convert_bits_to_random_int_below(uint32_t bits, uint32_t strict_max);

	int IRandom(int min, int max);
	//uint32_t MotherBits();

	// karl: no checks done
	unsigned int fast_random_int_below(unsigned int);

	inline void reseed(int seed1, int seed2);

//	void batch_BRandom(uint32_t* batch, int batch_size);
//
//	void nocheck_batch_BRandom(uint32_t* batch, int batch_size);

protected:
	uint32_t chunk_size() const { return SFMT_N * 4; }
	int left_in_chunk() const { return chunk_size() - ix; }

//	void Generate();
};

class RandomSFMT {
public:
	// BEWARE RandomSFMT(seed) and RandomSFMTY(seed,0) are not equivalent
	RandomSFMT(int seed = time(NULL)) : _impl(seed) {}
	RandomSFMT(int seed1, int seed2) : _impl(seed1, seed2) {}
	inline double random() { return _impl.Random(); }
	inline unsigned int random_int_below(unsigned int n) {
		return _impl.fast_random_int_below(n);
	}

	inline void reseed(int seed1, int seed2=0) {
		_impl.reseed(seed1, seed2);
	}


private:
	CRandomSFMT0_karl _impl;
};

//class RandomSFMT : public Random {
//public:
//	RandomSFMT(unsigned int seed = time(NULL)) : _impl(seed) {}
//	inline double random() { return _impl.Random(); }
//	inline int random_int_below(int n) {
//		return _impl.IRandom(0, n-1);
//	}
//
//private:
//	CRandomSFMT0_opt _impl;
//};

//// Subfunction for the sfmt algorithm
//static inline __m128i sfmt_recursion(__m128i const &a, __m128i const &b,
//__m128i const &c, __m128i const &d, __m128i const &mask) {
//    __m128i a1, b1, c1, d1, z1, z2;
//    b1 = _mm_srli_epi32(b, SFMT_SR1);
//    a1 = _mm_slli_si128(a, SFMT_SL2);
//    c1 = _mm_srli_si128(c, SFMT_SR2);
//    d1 = _mm_slli_epi32(d, SFMT_SL1);
//    b1 = _mm_and_si128(b1, mask);
//    z1 = _mm_xor_si128(a, a1);
//    z2 = _mm_xor_si128(b1, d1);
//    z1 = _mm_xor_si128(z1, c1);
//    z2 = _mm_xor_si128(z1, z2);
//    return z2;
//}

//inline void CRandomSFMT0_opt::Generate() {
//   // Fill state array with new random numbers
//   int i;
//   __m128i r, r1, r2;
//
//   r1 = state[SFMT_N - 2];
//   r2 = state[SFMT_N - 1];
//   for (i = 0; i < SFMT_N - SFMT_M; i++) {
//      r = sfmt_recursion(state[i], state[i + SFMT_M], r1, r2, mask);
//      state[i] = r;
//      r1 = r2;
//      r2 = r;
//   }
//   for (; i < SFMT_N; i++) {
//      r = sfmt_recursion(state[i], state[i + SFMT_M - SFMT_N], r1, r2, mask);
//      state[i] = r;
//      r1 = r2;
//      r2 = r;
//   }
//   ix = 0;
//}

//inline uint32_t CRandomSFMT0_opt::BRandom() {
//	// Output 32 random bits
//	uint32_t y;
//
//	if (ix >= SFMT_N * 4) {
//		Generate();
//	}
//	y = ((uint32_t*) state)[ix++];
//
//	if (UseMother)
//		y += MotherBits();
//	return y;
//}

inline CRandomSFMT0_karl::CRandomSFMT0_karl(int seed1, int seed2)
	: CRandomSFMT0(0)
{
	reseed(seed1, seed2);
}

inline void CRandomSFMT0_karl::reseed(int seed1, int seed2)
{
	int seeds[] = {seed1, seed2 };
	RandomInitByArray(seeds, 2);
}

inline uint32_t CRandomSFMT0_karl::BRandom() {
	// Output 32 random bits
	uint32_t y;

	if (ix >= SFMT_N * 4) {
		Generate();
	}
	y = ((uint32_t*) state)[ix++];
// karl
//	if (UseMother)
//		y += MotherBits();
	return y;
}

inline unsigned int CRandomSFMT0_karl::fast_random_int_below(unsigned int n)
{
//	// karl hack:
//	// Output random integer in the interval 0 <= x < n
//	// Slightly inaccurate if (n) is not a power of 2
//	// Assume 64 bit integers supported. Use multiply and shift method
//	uint32_t interval = (uint32_t)n; // Length of interval
//	uint64_t longran = (uint64_t) BRandom() * interval; // Random bits * interval
//	uint32_t iran = (uint32_t) (longran >> 32); // Longran / 2^32
//	return iran;
	return convert_bits_to_random_int_below(BRandom(), n);
}

inline unsigned int CRandomSFMT0_karl::convert_bits_to_random_int_below(uint32_t bits, uint32_t strict_max)
{
	uint64_t longran = (uint64_t)bits * strict_max;
	uint32_t iran = (uint32_t) (longran >> 32); // Longran / 2^32
	return iran;
}

//inline void CRandomSFMT0_karl::batch_BRandom(uint32_t* batch, int batch_size)
//{
//	const int CHUNK_SIZE = chunk_size();
//
//	int left = left_in_chunk();
//	if ( left > batch_size )
//		left = batch_size;
//	nocheck_batch_random_int_below(strict_max, batch, left);
//	batch += left;
//	batch_size -= left;
//
//	const int nb_chunks = batch_size / CHUNK_SIZE;
//	for (int i = 0; i < nb_chunks; ++i) {
//		Generate();
//		nocheck_batch_random_int_below(strict_max, batch, CHUNK_SIZE);
//		batch += CHUNK_SIZE;
//	}
//
//	int remaining = batch_size % CHUNK_SIZE;
//	if ( remaining > 0 ) {
//		Generate();
//		nocheck_batch_random_int_below(strict_max, batch, remaining);
//	}
//
//
//}
//
//inline void CRandomSFMT0_karl::nocheck_batch_BRandom(uint32_t* batch, int batch_size)
//{
//	for (int i = 0; i < batch_size; ++i)
//		batch[i] = ((uint32_t*) state)[ix++];
//}

inline int CRandomSFMT0_karl::IRandom(int min, int max) {
	// Output random integer in the interval min <= x <= max
	// Slightly inaccurate if (max-min+1) is not a power of 2
	if (max <= min) {
		if (max == min)
			return min;
		else
			return 0x80000000;
	}
	// Assume 64 bit integers supported. Use multiply and shift method
	uint32_t interval; // Length of interval
	uint64_t longran; // Random bits * interval
	uint32_t iran; // Longran / 2^32

	interval = (uint32_t) (max - min + 1);
	longran = (uint64_t) BRandom() * interval;
	iran = (uint32_t) (longran >> 32);
	// Convert back to signed and return result
	return (int32_t) iran + min;
}

//inline uint32_t CRandomSFMT0_opt::MotherBits() {
//   // Get random bits from Mother-Of-All generator
//   uint64_t sum;
//   sum =
//      (uint64_t)2111111111U * (uint64_t)MotherState[3] +
//      (uint64_t)1492 * (uint64_t)MotherState[2] +
//      (uint64_t)1776 * (uint64_t)MotherState[1] +
//      (uint64_t)5115 * (uint64_t)MotherState[0] +
//      (uint64_t)MotherState[4];
//   MotherState[3] = MotherState[2];
//   MotherState[2] = MotherState[1];
//   MotherState[1] = MotherState[0];
//   MotherState[4] = (uint32_t)(sum >> 32);       // Carry
//   MotherState[0] = (uint32_t)sum;               // Low 32 bits of sum
//   return MotherState[0];
//}

//inline int CRandomSFMT0_opt::IRandom(int min, int max) {
//	// Output random integer in the interval min <= x <= max
//	// Slightly inaccurate if (max-min+1) is not a power of 2
//	if (max <= min) {
//		if (max == min)
//			return min;
//		else
//			return 0x80000000;
//	}
//	// Assume 64 bit integers supported. Use multiply and shift method
//	uint32_t interval; // Length of interval
//	uint64_t longran; // Random bits * interval
//	uint32_t iran; // Longran / 2^32
//
//	interval = (uint32_t) (max - min + 1);
//	longran = (uint64_t) BRandom() * interval;
//	iran = (uint32_t) (longran >> 32);
//	// Convert back to signed and return result
//	return (int32_t) iran + min;
//}



#endif /* RANDOM_HPP_ */
