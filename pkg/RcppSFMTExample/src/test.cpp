#include <R.h>
#include <omp.h>

extern "C" {

#include "rngSFMTOpenMP.h"

void trySetRngParams(int *initType)
{
    setRngParams(*initType);
}

void tryGetStreams(int *pn, int *pm, double *out)
{
    int i, j, k;
    createRngArray(1, 0);
    k = 0;
    for (j=0; j < *pm; j++) {
        initializeRngInstance(0, j);
        for (i=0; i < *pn; i++) {
            out[k++] = getRngDouble(0);
        }
    }
    destroyRngArray();
}

void tryInterleavingSeq(int *pn, int *pm, int *ind, double *out)
{
    int i;
    createRngArray(*pn, 1);
    // store the numbers generated from different streams in the specified order
    for (i=0; i < *pm; i++) {
        out[i] = getRngDouble(ind[i]);
    }
    destroyRngArray();
}

void tryInterleavingThr(int *pn, int *pm, int *ind, double *out, int *nbThreads)
{
    int i;
    createRngArray(*pn, 1);
    omp_set_num_threads(*nbThreads);
    // store the numbers generated from different streams in the specified order
    #pragma omp parallel for schedule(dynamic)
    for (i=0; i < *pn; i++) {
        for (int j=0; j < *pm; j++) {
            if (ind[j] == i) {
                out[j] = getRngDouble(i);
            }
        }
    }
    destroyRngArray();
}

void trySumSequential(int *pn, int *length, double *out)
{
    int i, j;
    createRngArray(1, 0);
    // Compute the sum of the first length[i] numbers from stream i = 0, ..., n-1.
    for (i=0; i < *pn; i++) {
        initializeRngInstance(0, i);
        for (j=0; j < length[i]; j++) {
            out[i] += getRngDouble(0);
        }
    }
    destroyRngArray();
}

void trySumThreaded(int *pn, int *length, double *out, int *nbThreads)
{
    int i;
    createRngArray(*nbThreads, 0);
    omp_set_num_threads(*nbThreads);
    // Compute the sum of the first length[i] numbers from stream i = 0, ..., n-1.
    #pragma omp parallel for schedule(dynamic)
    for (i=0; i < *pn; i++) {
        int thread = omp_get_thread_num();
        initializeRngInstance(thread, i);
        for (int j=0; j < length[i]; j++) {
            out[i] += getRngDouble(thread);
        }
    }
    destroyRngArray();
}

} // extern "C"

