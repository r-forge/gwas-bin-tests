/*
 * rng_lecuyer.hpp
 *
 * An implementation of the RNG template interface for the rlecuyer implementation
 */
#ifndef _RCPPRANDOMSFMT_RNG_LECUYER_HPP_
#define _RCPPRANDOMSFMT_RNG_LECUYER_HPP_

extern "C" {
#include "RngStream.h"
}

class RNG_lecuyer {
public:

	RNG_lecuyer(int seed1, int seed2) {
		_stream = RngStream_CreateStream("RNG_lecuyer");
		set_seed(seed1, seed2);
	}
	~RNG_lecuyer() { RngStream_DeleteStream(&_stream);}

public: // ========== CLASS METHODS ==========
	double random() { return RngStream_RandU01(_stream); }
	int random_int_below(int n) {
		return RngStream_RandInt(_stream, 0, n);
	}
	void set_seed(int seed1, int seed2=0) {
		unsigned long seeds[6] = { seed1, seed2, seed1, seed2, seed1, seed2 };
		RngStream_SetSeed(_stream, seeds);
	}
private:
	RngStream _stream;

};

#endif /* _RCPPRANDOMSFMT_RNG_LECUYER_HPP_ */
