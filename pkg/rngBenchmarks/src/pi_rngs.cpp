#include <omp.h>
#include <rng_random_r.hpp>
#include <rng_sfmt.hpp>
#include "rng_lecuyer.hpp"
#include "rngOpenMP.h"
#include <Rcpp.h>

using namespace RcppRandomSFMT;
using namespace Rcpp;

template<class RNGT>
double pi_rng_threaded(unsigned long n, const int nb_threads, int seed )
{
	const int CHUNK_SIZE = 10000;

    RNGT* rngs[nb_threads];
    for (int i = 0; i < nb_threads; ++i)
    	rngs[i] = new RNGT(seed, i);

    omp_set_num_threads(nb_threads);

    const unsigned long  nb_chunks = n / (unsigned long)CHUNK_SIZE;

    double x = 0,y = 0;
    unsigned long inside = 0;

#pragma omp parallel for schedule(dynamic)
    for (unsigned long  i=0; i < nb_chunks; i++) {

    	int chunks_inside = 0;
    	int thread = omp_get_thread_num();
    	RNGT* rng = rngs[thread];
    	double x = 0,y = 0;
    	for (int j = 0; j < CHUNK_SIZE; j++) {
            x = rng->random()*2 - 1;
            y = rng->random()*2 - 1;
            if (x*x + y*y <= 1)
            	++chunks_inside;
    	}
    	// update the shared inside
#pragma omp atomic
    	inside += chunks_inside;

    }

    for (int i = 0; i < nb_threads; ++i)
    	delete rngs[i];

    double pih = (double)(inside)/(double)(n)*(double)(4);
    return pih;
}

double pi_rngOpenMP_threaded(unsigned long n, const int nb_threads, int seed )
{
	int CHUNK_SIZE = 10000;
    initializeRngArray(nb_threads);
    unsigned long inside = 0;

    omp_set_num_threads(nb_threads);

    unsigned long  nb_chunks = n / (unsigned long)CHUNK_SIZE;

    double x = 0,y = 0;
    //REprintf("n=%f\n", (float)n);
#pragma omp parallel for schedule(dynamic)
    for (unsigned long  i=0; i < nb_chunks; i++) {

    	int chunks_inside = 0;
    	int thread = omp_get_thread_num();
    	double x = 0,y = 0;
    	for (int j = 0; j < CHUNK_SIZE; j++) {
            x = getRngDouble(thread)*2 - 1;
            y = getRngDouble(thread)*2 - 1;
            if (x*x + y*y <= 1)
            	++chunks_inside;
    	}
    	// update the shared inside
#pragma omp atomic
    	inside += chunks_inside;

    }
    destroyRngArray();
    return (double)(inside)/(double)(n)*(double)(4);
}

RcppExport SEXP _pi_rngOpenMP_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed) {
	double res = pi_rngOpenMP_threaded( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_random_r_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed) {
	double res = pi_rng_threaded<RNGrandomR>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_sfmt_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed) {
	double res = pi_rng_threaded<RNGSFMT>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_lecuyer_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed) {
	double res = pi_rng_threaded<RNG_lecuyer>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed) );
	return wrap(res);
}
