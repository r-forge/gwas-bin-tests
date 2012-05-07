#include <R_ext/Rdynload.h>
#include "rngSFMTOpenMP.h"

void R_init_rngVerify(DllInfo *info)
{
    createRngArray = (void (*) (int, int)) R_GetCCallable("rngSFMTOpenMP", "createRngArray");
    initializeRngInstance = (void (*) ()) R_GetCCallable("rngSFMTOpenMP", "initializeRngInstance");
    destroyRngArray = (void (*) ()) R_GetCCallable("rngSFMTOpenMP", "destroyRngArray");
    getRngDouble = (double (*) (int)) R_GetCCallable("rngSFMTOpenMP", "getRngDouble");
    setRngParams = (void (*) (int)) R_GetCCallable("rngSFMTOpenMP", "setRngParams");
}

