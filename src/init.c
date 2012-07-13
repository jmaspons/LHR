/*
 *  Native routines registration, as per "Writing R extensions".
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static const R_ExternalMethodDef ExternalEntries[] = {
    {"actuar_do_random", (DL_FUNC) &actuar_do_random, -1},
    {"actuar_do_dpq", (DL_FUNC) &actuar_do_dpq, -1},
    {NULL, NULL, 0}
};

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
}
