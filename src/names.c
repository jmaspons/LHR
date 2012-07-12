/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Table of functions internal to the package. First element is an
 *  argument to one of actuar_do_dpq or actuar_do_random, functions
 *  callable from .External(); second element is the C function
 *  actually called; third element is a code used in the latter.
 *
 *  Idea taken from R sources (see .../src/main/names.c).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <Rinternals.h>
#include "actuar.h"

FUNTAB fun_tab[] = {
    /* DENSITY, CUMULATIVE PROBABILITY AND QUANTILE FUNCTIONS,
     * RAW AND LIMITED MOMENTS */
    /* Three parameter distributions */
    {"dbetabinom",       actuar_do_dpq3,         1},
    {"pbetabinom",       actuar_do_dpq3,         2},
    {"qbetabinom",       actuar_do_dpq3,         3},
    {"mbetabinom",       actuar_do_dpq3,         4},
    {"dbetanbinom",      actuar_do_dpq3,         5},
    {"pbetanbinom",      actuar_do_dpq3,         6},
    {"qbetanbinom",      actuar_do_dpq3,         7},
    {"mbetanbinom",      actuar_do_dpq3,         8},

    /* RANDOM NUMBERS FUNCTIONS */
    /* Three parameter distributions */
    {"rbetabinom",       actuar_do_random3,     1},
    {"rbetanbinom",      actuar_do_random3,     2},
    {0, 0, 0}
};
