/*
 * rng_rand.hpp
 *
 * An implementation of the RNG template interface using the system rand() function
 * This is mostly for benchmarking reasons because the system rand() is usually neither fast nor good.
 * N.B: it is not reentrant so it can not be used in multithread programs
 *
 */
#ifndef _RCPPRANDOMSFMT_RNG_RAND_HPP_
#define _RCPPRANDOMSFMT_RNG_RAND_HPP_

namespace RcppRandomSFMT {

#include <stdlib.h>
#include "rng.hpp"

class RNGrand {
public:
	RNGrand(int seed=time(NULL), int seed2 = 0) { srand(seed); }
	~RNGrand() {}

public: // ========== CLASS METHODS ==========
	double random() {  rand()/(RAND_MAX+1.0);  }
	int random_int_below(int n) { return (int) ( ((double)n)* ( rand()/(RAND_MAX+1.0) ) );}
};

} // end of namespace

#endif /* _RCPPRANDOMSFMT_RNG_RAND_HPP_ */
