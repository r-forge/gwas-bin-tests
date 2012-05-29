#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;
#include "SFMT.hpp"


// Rcpp_MODULE does not like default args
// so we work around
void _sfmt_seet_seeds_from_R(RNG::SFMT* self) {
	self->set_seeds_from_R();
}

// ===== CLASS METHODS ====
// a priori, Rcpp_MODULE does not support class methods
// so we work-around this as well

int _GET_MERSENNE_EXPONENT(RNG::SFMT* self) {
	return RNG::SFMT::GET_MERSENNE_EXPONENT();
}

string _ID(RNG::SFMT* self) {
	return RNG::SFMT::ID();
}

// ==== useful instance methods ====

double _nextUnsignedRandomInteger(RNG::SFMT* self) {
	return double(self->nextRandomInteger());
}

IntegerVector _randomIntegers(RNG::SFMT* self, int n) {
	IntegerVector v(n);
	for (int i = 0; i < n; ++i)
		v[i] = self->nextRandomInteger();
	return v;
}

NumericVector _unsignedRandomIntegers(RNG::SFMT* self, int n) {
	NumericVector v(n);
	for (int i = 0; i < n; ++i)
		v[i] = self->nextRandomInteger();
	return v;
}

NumericVector _randomDoubles(RNG::SFMT* self, int n) {
	NumericVector v(n);
	for (int i = 0; i < n; ++i)
		v[i] = self->nextRandomDouble();
	return v;
}

void _set_seed_array(RNG::SFMT* self, vector<uint32_t> seeds) {
	self->set_seed_array(&seeds[0], seeds.size());
}

RCPP_MODULE(sfmt){
	using namespace Rcpp;
	using namespace RNG;

//	function("SFMT_MERSENNE_EXPONENT",
//			&RNG::SFMT::GET_MERSENNE_EXPONENT, "Get the Mersenne Exponent used by this SFMT implementation");
//
//	function("SFMT_ID", &RNG::SFMT::ID, "Get the SFMT identification string for this implementation");

    // we expose the class std::vector<double> as "vec" on the R side
    class_<SFMT>("SFMT")

    // exposing the constructors
    .constructor("initialize the RNG using the R current seed, cf set_seeds_from_R()")
    .constructor<int>("initialize the RNG using the seed given as argument")
    .constructor<int, int>("initialize the RNG using the pair of seeds given as argument")

    // exposing member functions
    .method( "nextUnsignedRandomInteger", &_nextUnsignedRandomInteger, "generate an usigned random 32 bits integer as a double")
    .method( "nextRandomInteger", &SFMT::nextRandomInteger, "generate a random 32 bits integer")
    .method( "nextRandomDouble", &SFMT::nextRandomDouble, "generate a random double in [0,1[")
    .method( "set_seeds_from_R", &_sfmt_seet_seeds_from_R, "set the seed of the RNG to the current R seed (cf set.seed() )")
    .method( "set_seed", &SFMT::set_seed, "(re)set the seed of the RNG to the supplied seed")
    .method( "set_seeds", &SFMT::set_seeds, "(re)set the seed of the RNG to the supplied pair of seeds")
    .method( "set_seed_array", &_set_seed_array, "(re)set the seed of the RNG to the supplied vector of seeds")

    .method("unsignedRandomIntegers", &_unsignedRandomIntegers, "generate a vector of random unsigned 32 bits integers as a numeric vector")
    .method("randomIntegers", &_randomIntegers, "generate a vector of random 32 bits integers")
    .method("randomDoubles", &_randomDoubles, "generate a vector of random doubles")

    .method("get_mersenne_exponent", &_GET_MERSENNE_EXPONENT)
    .method("id", &_ID)

    ;
}
