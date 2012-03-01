/*
 * rng_sfmt.hpp
 *
 * An implementation of the RNG template interface using modified Randomc CrandomSFMT0 Generator
 */
#ifndef _RCPPRANDOMSFMT_RNG_SFMT_HPP_
#define _RCPPRANDOMSFMT_RNG_SFMT_HPP_

#include "rng.hpp"
#include "crandomsfmt0_modified.hpp"

namespace RcppRandomSFMT {

class RNGSFMT {
public:
	RNGSFMT(int seed = time(NULL)) : _impl(seed) {}
	RNGSFMT(int seed1, int seed2) : _impl(seed1, seed2) {}
	~RNGSFMT() {}

public: // ========== CLASS METHODS ==========
	double random() { return _impl.Random(); }
	int random_int_below(int n) {
			return _impl.fast_random_int_below(n);
	}
	void reseed(int seed1, int seed2=0) {
			_impl.reseed(seed1, seed2);
	}
private:
	CRandomSFMT0_modified _impl;

};

} // end of namespace

#endif /* _RCPPRANDOMSFMT_RNG_SFMT_HPP_ */
