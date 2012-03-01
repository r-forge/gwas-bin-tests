/*
 * crandomsfmt0_modified.hpp
 *
 *  a slightly modified version of the CRandomSFMT0 generator of the Randomc library
 *
 *  It was modified to further increase speed on my platform, by inlining some
 *  methods and removing the feature to use the Mother Generator
 *
 *  Author: karl.forner@gmail.com
 *
 */

#ifndef _RCPPRANDOMSFMT_CRANDOMSFMT0_MODIFIED_HPP_
#define _RCPPRANDOMSFMT_CRANDOMSFMT0_MODIFIED_HPP_

#include <sfmt.h>
#include <sfmt_inline.hpp>

class CRandomSFMT0_modified : public CRandomSFMT0 {
public: // ===== LIFECYCLE =====
	CRandomSFMT0_modified(int seed) : CRandomSFMT0(seed) {}
	CRandomSFMT0_modified(int seed1, int seed2) : CRandomSFMT0(0) { reseed(seed1, seed2); }

public: // ===== OVERRIDDEN  CRandomSFMT0 METHODS =====
	// Give 32 random bits
	// inlined the original CRandomSFMT::BRandom method in sfmt.cpp
	// and removed the if (UseMother) test
	uint32_t BRandom();

	// give random number between 0 and 1
	// modified to remove the if (UseMother) test
	 double Random();

	// inlined the unmodified original CRandomSFMT::IRandom method in sfmt.cpp
	int IRandom(int min, int max);

public: // ===== NEW PUBLIC INSTANCE METHODS =====
	// transform a random number (32 random bits) into a random number < strict_max
	// I do not remember where did I take this from
	unsigned int convert_bits_to_random_int_below(uint32_t bits, uint32_t strict_max)
	{
		uint64_t longran = (uint64_t)bits * strict_max;
		uint32_t iran = (uint32_t) (longran >> 32); // Longran / 2^32
		return iran;
	}


	// fast way to use the generator to output a random number x such as: 0 <= x < max
	// i.e. it is a faster version of IRandom(min, max) when min==0
	// and
	unsigned int fast_random_int_below(unsigned int max) {
		return convert_bits_to_random_int_below(BRandom(), max);
	}

	// reseed the generator with a pair of seeds
	inline void reseed(int seed1, int seed2);

};


inline void CRandomSFMT0_modified::reseed(int seed1, int seed2)
{
	int seeds[] = {seed1, seed2 };
	RandomInitByArray(seeds, 2);
}

inline double CRandomSFMT0_modified::Random() {
   // Output random floating point number
   if (ix >= SFMT_N*4-1) {
      // Make sure we have at least two 32-bit numbers
      Generate();
   }
   uint64_t r = *(uint64_t*)((uint32_t*)state+ix);
   ix += 2;
//   if (UseMother) {
//      // We need 53 bits from Mother-Of-All generator
//      // Use the regular 32 bits and the the carry bits rotated
//      uint64_t r2 = (uint64_t)MotherBits() << 32;
//      r2 |= (MotherState[4] << 16) | (MotherState[4] >> 16);
//      r += r2;
//   }
   // 53 bits resolution:
   // return (int64_t)(r >> 11) * (1./(67108864.0*134217728.0)); // (r >> 11)*2^(-53)
   // 52 bits resolution for compatibility with assembly version:
   return (int64_t)(r >> 12) * (1./(67108864.0*67108864.0));  // (r >> 12)*2^(-52)
}

inline uint32_t CRandomSFMT0_modified::BRandom() {
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


// unmodified code from sfmt.cpp
inline int CRandomSFMT0_modified::IRandom(int min, int max) {
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



#endif /* _RCPPRANDOMSFMT_CRANDOMSFMT0_MODIFIED._HPP_ */
