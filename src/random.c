/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to generate variates of some probability laws not in
 *  base R. Function .External() calls actuar_do_random() with
 *  arguments:
 *
 *      1. the name of the distribution from which to simulate, with
 *         an "r" prepended to it (e.g. "rpareto");
 *      2. the number of variates;
 *    3:x. the parameters of the distribution.
 *
 *  Function actuar_do_random() will extract the name of the
 *  distribution, look up in table fun_tab defined in names.c which of
 *  actuar_do_random{1,2,3,4} should take care of the simulation and
 *  dispatch to this function. In turn, functions
 *  actuar_do_random{1,2,3,4} call function rdist() to get actual
 *  variates from distribution "dist".
 *
 *  This scheme is essentially what is used in base R (see files
 *  .../src/main/random.c, .../src/main/names.c in R sources).
 *
 * To add a new distribution: write an rdist() function, add an entry
 * in names.c and in the definition of the corresponding
 * actuar_do_random{1,2,3,4} function, declare the function in
 * actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much indirect help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"

/* Functions for three parameter distributions */
static Rboolean random3(double (*f) (), double *a, int na,
                        double *b, int nb, double *c, int nc,
                        double *x, int n)
{
    double ai, bi, ci;
    int i;
    Rboolean naflag = FALSE;
    for (i = 0; i < n; i++)
    {
        ai = a[i % na];
        bi = b[i % nb];
        ci = c[i % nc];
        x[i] = f(ai, bi, ci);
        if (ISNAN(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RAND3(num, fun) \
        case num: \
            naflag = random3(fun, REAL(a), na, REAL(b), nb, REAL(c), nc, REAL(x), n); \
            break


SEXP actuar_do_random3(int code, SEXP args)
{
    SEXP x, a, b, c;
    int i, n, na, nb, nc;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
            !isNumeric(CADR(args)) ||
            !isNumeric(CADDR(args)) ||
            !isNumeric(CADDDR(args)))
        error(_("invalid arguments"));

    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1)
    {
        n = asInteger(CAR(args));
        if (n == NA_INTEGER || n < 0)
            error(_("invalid arguments"));
    }
    else
        n = LENGTH(CAR(args));

    /* If n == 0, return numeric(0) */
    PROTECT(x = allocVector(REALSXP, n));
    if (n == 0)
    {
        UNPROTECT(1);
        return(x);
    }

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    nc = LENGTH(CADDDR(args));
    if (na < 1 || nb < 1 || nc < 1)
    {
        for (i = 0; i < n; i++)
            REAL(x)[i] = NA_REAL;
        warning(_("NAs produced"));
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
        Rboolean naflag = FALSE;
        PROTECT(a = coerceVector(CADR(args), REALSXP));
        PROTECT(b = coerceVector(CADDR(args), REALSXP));
        PROTECT(c = coerceVector(CADDDR(args), REALSXP));
        GetRNGstate();
        switch (code)
        {
            RAND3(1, rbetabinom);
            RAND3(2, rbetanbinom);
        default:
            error(_("internal error in actuar_do_random3"));
        }
        if (naflag)
            warning(R_MSG_NA);
        PutRNGstate();
        UNPROTECT(3);
    }
    UNPROTECT(1);
    return x;
}

/* Main function, the only one used by .External(). */
SEXP actuar_do_random(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to actuar_do_random{1,2,3,4} */
    for (i = 0; fun_tab[i].name; i++)
    {
        if (!strcmp(fun_tab[i].name, name))
            return fun_tab[i].cfun(fun_tab[i].code, CDR(args));
    }

    /* No dispatch is an error */
    error(_("internal error in actuar_do_random"));

    return args;                /* never used; to keep -Wall happy */
}
