/*
 * rng_mrg32k5a.hpp
 *
 * An implementation of the RNG template interface using the rngOpenMP implementation of mrg32k5a
 */
#ifndef _RCPPRANDOMSFMT_RNG_MRG32K5A_HPP_
#define _RCPPRANDOMSFMT_RNG_MRG32K5A_HPP_

#include "classRandom.h"

class RNG_MRG32K5A {
public:

	RNG_MRG32K5A(int seed1, int seed2 = 0): _impl()  {
		set_seed(seed1, seed2);
	}
	~RNG_MRG32K5A() {}

public: // ========== CLASS METHODS ==========
	double random() { return _impl.getDouble(); }
	int random_int_below(int n) {
		return int(random()*n);
	}
	void set_seed(int seed1, int seed2=0) {
		double state[10];
		state[0] = seed1;
		state[1] = seed2;
		_impl.putState(&(state[0]));
	}
private:
	ClassRandom _impl;

};

#endif /* _RCPPRANDOMSFMT_RNG_LECUYER_HPP_ */
