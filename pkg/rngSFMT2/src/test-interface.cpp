#include <R.h>

#include "SFMT_include_all.hpp"
//#include "Sfmt_607_NORMAL.hpp"
//#include "Sfmt_19937_NORMAL.hpp"
//RNGSFMT *rng;

extern "C" {

void sfmt_607_test(int *seed, int *nb, double *res)
{
	SFMT_607_NORMAL::RNG s(*seed);
	double sum = 0;
	int n = *nb;
	for (int i = 0; i < n; ++i) {
		sum += s.rand32();
	}
	*res = sum;
}


void sfmt_19937NORMAL(int *seed, int *nb, double *res)
{
	SFMT_19937_NORMAL::RNG s(*seed);
	double sum = 0;
	int n = *nb;
	for (int i = 0; i < n; ++i) {
		sum += s.rand32();
	}
	*res = sum;
}

void sfmt_19937SSE2(int *seed, int *nb, double *res)
{
	SFMT_19937_SSE2::RNG s(*seed);
	double sum = 0;
	int n = *nb;
	for (int i = 0; i < n; ++i) {
		sum += s.rand32();
	}
	*res = sum;
}

} // extern "C"

