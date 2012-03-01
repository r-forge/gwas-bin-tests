#include <R.h>
#include "rngOpenMP.h"

void testInterleaving(int *pn, int *pm, int *ind, double *out)
{
    int i;
    initializeRngArray(*pn);
    // store the numbers generated from different streams in the specified order
    for (i=0; i < *pm; i++) {
        out[i] = getRngDouble(ind[i]);
    }
    destroyRngArray();
}

void testSequential(int *pn, int *length, double *out)
{
    int n, i, j;
    n = *pn;
    initializeRngArray(n);
    // Compute the sum of the first length[i] numbers from stream i, i = 0, ..., n-1.
    for (i=0; i < n; i++) {
        for (j=0; j < length[i]; j++) {
            out[i] += getRngDouble(i);
        }
    }
    destroyRngArray();
}

void testThreaded(int *pn, int *length, double *out)
{
    int n, i, j;
    n = *pn;
    initializeRngArray(n);
    // Compute the sum of the first length[i] numbers from stream i, i = 0, ..., n-1.
    #pragma omp parallel for private(j)
    for (i=0; i < n; i++) {
        for (j=0; j < length[i]; j++) {
            out[i] += getRngDouble(i);
        }
    }
    destroyRngArray();
}

