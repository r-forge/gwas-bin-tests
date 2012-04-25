#include <R.h>

#include "SFMT.hpp"
//#include "Sfmt_607_NORMAL.hpp"
//#include "Sfmt_19937_NORMAL.hpp"
//RNGSFMT *rng;

extern "C" {

void sfmt_test(int *seed, int *nb, double *res)
{
	RNG::SFMT s(*seed);
	double sum = 0;
	int n = *nb;
	for (int i = 0; i < n; ++i) {
		sum += s.rand32();
	}
	*res = sum;
}


} // extern "C"

