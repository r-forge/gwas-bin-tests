#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;
#include "SFMT.hpp"

//9937_NORMAL.hpp"
//RNGSFMT *rng;
extern "C" {


RcppExport SEXP sfmt_test(int *seed, int *nb, double *res)
{
	RNG::SFMT s(*seed);
	double sum = 0;
	int n = *nb;
	for (int i = 0; i < n; ++i) {
		sum += s.rand32();
	}
	*res = sum;
}

RcppExport SEXP W_fill_array32(SEXP _nb, SEXP _seed) {
	int seed = as<int>(_seed);
	int nb = as<int>(_nb);
	vector<uint32_t> v;
	v.resize(nb);
	RNG::SFMT s(seed);
	s.rand32(v);

	vector<double> res;
	res.resize(nb);
	for (int i = 0; i < nb; ++i)
		res[i] = double(v[i]);


	return wrap(res);
}


RcppExport SEXP  MERSENNE_EXPONENT() {
	return wrap(RNG::SFMT::GET_MERSENNE_EXPONENT());
}


RcppExport SEXP  sfmt_printid() {
	Rprintf("SFMT ID=%s\n", RNG::SFMT::ID().c_str());
	return R_NilValue;
}


} // extern "C"

