#include <R.h>
#include <omp.h>
#include "rngSFMTOpenMP.h"

void pi_rng_threaded(int *pn, double *out, int *nbThreads)
{
    const int CHUNK_SIZE = 10000;
    int i, j, thread, chunks_inside;
    double x, y;
    unsigned long inside = 0;

    omp_set_num_threads(*nbThreads);
    createRngArray(*nbThreads, 0);

    const unsigned long  nb_chunks = *pn / (unsigned long)CHUNK_SIZE;


    #pragma omp parallel for private(thread, j, x, y, chunks_inside)
    for (i=0; i < nb_chunks; i++) {
        thread = omp_get_thread_num();
        initializeRngInstance(thread, i);
        chunks_inside = 0;
        for (j = 0; j < CHUNK_SIZE; j++) {
            x = getRngDouble(thread);
            y = getRngDouble(thread);
            if (x*x + y*y <= 1)
                ++chunks_inside;
        }
        // update the shared inside
        #pragma omp atomic
        inside += chunks_inside;
    }
    destroyRngArray();

    *out = (double)(inside)/(double)(nb_chunks * CHUNK_SIZE)*4.0;
}

