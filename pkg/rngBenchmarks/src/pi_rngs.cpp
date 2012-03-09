#include <omp.h>
#include <rng_random_r.hpp>
#include <rng_sfmt.hpp>
#include "rng_lecuyer.hpp"
//#include "rngOpenMP.h"
#include "rng_mrg32k5a.hpp"
#include <Rcpp.h>

using namespace RcppRandomSFMT;
using namespace Rcpp;

template<class RNGT> int process_chunk(RNGT* rng, const int nb) {
	int inside = 0;
	double x = 0,y = 0;
	for (int j = 0; j < nb; j++) {
        x = rng->random();
        y = rng->random();
        if (x*x + y*y <= 1)
        	++inside;
	}
	return inside;
}

template<class RNGT>
double pi_rng_threaded(unsigned long n, const int nb_threads, int seed, int chunk_size )
{
    RNGT* rngs[nb_threads];
    for (int i = 0; i < nb_threads; ++i)
    	rngs[i] = new RNGT(seed, i);

    omp_set_num_threads(nb_threads);

    const unsigned long nb_chunks = n / (unsigned long)chunk_size;
    const int remaining = n % (unsigned long)chunk_size;

    unsigned long inside = 0;

#pragma omp parallel for schedule(dynamic)
    for (unsigned long  i=0; i < nb_chunks; i++) {
    	RNGT* rng = rngs[ omp_get_thread_num() ];
    	rng->set_seed( seed, i );
    	int chunks_inside = process_chunk<RNGT>(rng, chunk_size);
    	// update the shared inside
#pragma omp atomic
    	inside += chunks_inside;
    }
    // compute the left-over in the main thread, no need to synchronyze
    if ( remaining > 0 ) {
    	RNGT* rng = rngs[ 0 ];
    	rng->set_seed( seed, nb_chunks );
    	inside += process_chunk<RNGT>(rng, remaining);
    }

    for (int i = 0; i < nb_threads; ++i)
    	delete rngs[i];

    return (double)(inside)/(double)(n)*(double)(4);
}
/*
double pi_rngOpenMP_threaded(unsigned long n, const int nb_threads, int seed, , int chunk_size  )
{
	int CHUNK_SIZE = 10000;
    initializeRngArray(nb_threads);
    unsigned long inside = 0;

    omp_set_num_threads(nb_threads);

    unsigned long  nb_chunks = n / (unsigned long)chunk_size;

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
*/

//RcppExport SEXP _pi_rngOpenMP_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
//	double res = pi_rngOpenMP_threaded( as<unsigned long>(_n)
//			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
//	return wrap(res);
//}

RcppExport SEXP _pi_rng_random_r_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
	double res = pi_rng_threaded<RNGrandomR>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_sfmt_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
	double res = pi_rng_threaded<RNGSFMT>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_lecuyer_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
	double res = pi_rng_threaded<RNG_lecuyer>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
	return wrap(res);
}

RcppExport SEXP _pi_rng_mrg32k5a_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
	double res = pi_rng_threaded<RNG_MRG32K5A>( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
	return wrap(res);
}
