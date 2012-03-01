/*
 * rng.hpp
 *
 * A base class for the Random generators provided
 *
 * N.B: This class does not provide virtual methods, so you can not
 * use polymorphism on it, for efficiency reasons.
 *
 * Its only purpose is to describe the common interface used by derived classes
 * so that the derived classes may be used in templates for optimal efficiency
 *
 */
#ifndef _RCPPRANDOMSFMT_RNG_HPP_
#define _RCPPRANDOMSFMT_RNG_HPP_

namespace RcppRandomSFMT {

#include <time.h>

class RNG {
public:
	RNG(int seed=time(NULL), int seed2 = 0) {  }
	~RNG() {}

public: // ========== CLASS METHODS ==========
	double random() { return 0.0; }
	int random_int_below(int n) { return 0;}
};

} // end of namespace

#endif /* _RCPPRANDOMSFMT_RNG_HPP_ */
