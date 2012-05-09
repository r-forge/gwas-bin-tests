#include <omp.h>
#include <Rcpp.h>
#include <SFMT.hpp>

using namespace RNG;
using namespace Rcpp;

int process_chunk(SFMT* rng, const int nb) {
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


double pi_rng_threaded(unsigned long n, const int nb_threads, int seed, int chunk_size )
{
	SFMT* rngs[nb_threads];
    for (int i = 0; i < nb_threads; ++i)
    	rngs[i] = new SFMT(seed, i);

    omp_set_num_threads(nb_threads);

    const unsigned long nb_chunks = n / (unsigned long)chunk_size;
    const int remaining = n % (unsigned long)chunk_size;

    unsigned long inside = 0;

#pragma omp parallel for schedule(dynamic)
    for (unsigned long  i=0; i < nb_chunks; i++) {
    	SFMT* rng = rngs[ omp_get_thread_num() ];
    	rng->set_seeds( seed, i );
    	int chunks_inside = process_chunk(rng, chunk_size);
    	// update the shared inside
#pragma omp atomic
    	inside += chunks_inside;
    }
    // compute the left-over in the main thread, no need to synchronyze
    if ( remaining > 0 ) {
    	SFMT* rng = rngs[ 0 ];
    	rng->set_seeds( seed, nb_chunks );
    	inside += process_chunk(rng, remaining);
    }

    for (int i = 0; i < nb_threads; ++i)
    	delete rngs[i];

    return (double)(inside)/(double)(n)*(double)(4);
}

RcppExport SEXP _pi_rng_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
	double res = pi_rng_threaded( as<unsigned long>(_n)
			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
	return wrap(res);
}

//RcppExport SEXP _pi_rngOpenMP_threaded(SEXP _n, SEXP _nb_threads, SEXP _seed, SEXP _chunk_size) {
//	double res = pi_rngOpenMP_threaded( as<unsigned long>(_n)
//			, as<int>(_nb_threads), as<int>(_seed), as<int>(_chunk_size) );
//	return wrap(res);
//}
