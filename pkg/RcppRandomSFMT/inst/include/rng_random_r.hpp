/*
 * rng_radom_r.hpp
 *
 * An implementation of the RNG template interface using the stdlib random_r
 * reentrant RNG
 */
#ifndef _RCPPRANDOMSFMT_RNG_RANDOM_R_HPP_
#define _RCPPRANDOMSFMT_RNG_RANDOM_R_HPP_

#include <cstdlib>
#include <time.h>
#include "rng.hpp"

namespace RcppRandomSFMT {


const int STATE_LENGTH = 256;

class RNGrandomR {
public:
	RNGrandomR(int seed=time(0), int seed2 = 0) {
		set_seed(seed);
	}
	~RNGrandomR() {}

public: // ========== CLASS METHODS ==========

	// N.B: in any case seed2 is IGNORED
	void set_seed(int seed1, int seed2=0) {
		_random_data.state = 0;
		initstate_r(seed1, _statebuf, STATE_LENGTH, &_random_data);
	}

	double random() {
		int32_t res;
		random_r(&_random_data, &res);
		return res / (RAND_MAX+1.0);
	}
	int random_int_below(int n) {
		return (int) ( ((double)n)* ( this->random()/(RAND_MAX+1.0) ) );
	}

private:
	struct random_data _random_data;
	char _statebuf[STATE_LENGTH];
};

}// end of namespace

#endif /* _RCPPRANDOMSFMT_RNG_R_HPP_ */
