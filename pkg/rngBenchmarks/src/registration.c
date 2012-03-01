#include <R_ext/Rdynload.h>
#include "rngOpenMP.h"

void R_init_rngBenchmarks(DllInfo *info)
{
    initializeRngArray = (void (*) (int)) R_GetCCallable("rngOpenMP", "initializeRngArray");
    getRngDouble = (double (*) (int)) R_GetCCallable("rngOpenMP", "getRngDouble");
    destroyRngArray = (void (*) ()) R_GetCCallable("rngOpenMP", "destroyRngArray");
}


