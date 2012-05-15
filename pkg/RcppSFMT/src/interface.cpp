#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;
#include "SFMT.hpp"

//9937_NORMAL.hpp"
//RNGSFMT *rng;
extern "C" {


//RcppExport SEXP sfmt_test(int *seed, int *nb, double *res)
//{
//	RNG::SFMT s(*seed);
//	double sum = 0;
//	int n = *nb;
//	for (int i = 0; i < n; ++i) {
//		sum += s.rand32();
//	}
//	*res = sum;
//}

//RcppExport SEXP W_fill_array32(SEXP _nb, SEXP _seed) {
//	int seed = as<int>(_seed);
//	int nb = as<int>(_nb);
//
//	RNG::SFMT s(seed);
//	vector<double> res;
//	res.resize(nb);
//	for (int i = 0; i < nb; ++i)
//		res[i] = double(s.nextRandomInteger());
//
//	return wrap(res);
//}

RcppExport SEXP W_fill_array32_seeds(SEXP _nb, SEXP _seeds) {
	IntegerVector seeds1(_seeds);

	vector<uint32_t> seeds(seeds1.begin(), seeds1.end());
	int nb = as<int>(_nb);

	RNG::SFMT s;
	if ( seeds.size() == 1)
		s.set_seed(seeds[0]);
	else
		s.set_seeds(&seeds[0], seeds.size());

	vector<double> res;
	res.resize(nb);
	for (int i = 0; i < nb; ++i)
		res[i] = double(s.nextRandomInteger());

	return wrap(res);
}


RcppExport SEXP  MERSENNE_EXPONENT() {
	return wrap(RNG::SFMT::GET_MERSENNE_EXPONENT());
}


RcppExport SEXP  SFMT_ID() {
	return wrap( RNG::SFMT::ID() );
}


} // extern "C"

