/*  ===== FROM: actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability density, cumulative probability
 *  quantile functions and moment generating functions, raw moments
 *  and limited moments for some probability laws not in base R (or
 *  those quantities not provided in base R). Function .External()
 *  calls actuar_do_dpq() with arguments:
 *
 *       1. the name of the distribution, with a "d", a "p" or "q"
 *          prepended to it (e.g. "dpareto", "pburr");
 *       2. the value(s) where the function is to be evaluated;
 *     3:x. the parameters of the distribution (including the order of
 *          the limited moment for lev*);
 *     x+1. whether to return the lower or upper tail probability or
 *          quantile (p* and q* only); see note below for m* and lev*
 *          functions;
 *     x+2. whether to return probability in log scale or the cumulant
 *          generating function (d*, p*, q* and mgf* only).
 *
 *  Function actuar_do_dpq() will extract the name of the distribution, look
 *  up in table fun_tab defined in names.c which of actuar_do_dpq{1,2,3,4}
 *  should take care of the calculation and dispatch to this function.
 *  In turn, functions actuar_do_dpq{1,2,3,4} call function
 *  {d,p,q,m,lev,mgf}dist() to get actual values from distribution
 *  "dist".
 *
 *  Note: the m* and lev* functions came later in the process. In
 *  order to easily fit them into this system, I have decided to leave
 *  an unused 'give_log' argument in the C definitions of these
 *  functions. Otherwise, this would have required defining functions
 *  dpq{1,2,3,4,5}_0() below.
 *
 *  Functions therein are essentially identical to those found in
 *  .../src/main/arithmetic.c of R sources with a different naming
 *  scheme.
 *
 *  To add a new distribution: write a {d,p,q,m,lev,mgf}dist()
 *  function, add an entry in names.c and in the definition of the
 *  corresponding actuar_do_dpq{1,2,3,4} function, declare the function in
 *  actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much indirect help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"


/* Functions for three parameter distributions */
#define if_NA_dpq3_set(y, x, a, b, c)                                       \
        if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c)) y = NA_REAL; \
        else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c)) y = R_NaN;

#define mod_iterate3(n1, n2, n3, n4, i1, i2, i3, i4)    \
        for (i = i1 = i2 = i3 = i4 = 0; i < n;          \
             i1 = (++i1 == n1) ? 0 : i1,                \
             i2 = (++i2 == n2) ? 0 : i2,                \
             i3 = (++i3 == n3) ? 0 : i3,                \
             i4 = (++i4 == n4) ? 0 : i4,                \
             ++i)

static SEXP dpq3_1(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, n, nx, na, nb, nc,
        sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb), sco = OBJECT(sc);
    double xi, ai, bi, ci, *x, *a, *b, *c, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ3                                              \
    if (!isNumeric(sx) || !isNumeric(sa) ||                     \
        !isNumeric(sb) || !isNumeric(sc))                       \
        error(_("invalid arguments"));                          \
                                                                \
    nx = LENGTH(sx);                                            \
    na = LENGTH(sa);                                            \
    nb = LENGTH(sb);                                            \
    nc = LENGTH(sc);                                            \
    if ((nx == 0) || (na == 0) || (nb == 0) || (nc == 0))       \
        return(allocVector(REALSXP, 0));                        \
    n = nx;                                                     \
    if (n < na) n = na;                                         \
    if (n < nb) n = nb;                                         \
    if (n < nc) n = nc;                                         \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sc = coerceVector(sc, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
    c = REAL(sc);                                               \
    y = REAL(sy)

    SETUP_DPQ3;

    i_1 = asInteger(sI);

    mod_iterate3(nx, na, nb, nc, ix, ia, ib, ic)
    {
        xi = x[ix];
        ai = a[ia];
        bi = b[ib];
        ci = c[ic];
        if_NA_dpq3_set(y[i], xi, ai, bi, ci)
        else
        {
            y[i] = f(xi, ai, bi, ci, i_1);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

#define FINISH_DPQ3                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
    else if (n == nb) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));  \
        SET_OBJECT(sy, sbo);                    \
    }                                           \
    else if (n == nc) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sc)));  \
        SET_OBJECT(sy, sco);                    \
    }                                           \
    UNPROTECT(5)

    FINISH_DPQ3;

    return sy;
}

static SEXP dpq3_2(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, n, nx, na, nb, nc,
        sxo = OBJECT(sx), sao = OBJECT(sa),
        sbo = OBJECT(sb), sco = OBJECT(sc);
    double xi, ai, bi, ci, *x, *a, *b, *c, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQ3;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);

    mod_iterate3(nx, na, nb, nc, ix, ia, ib, ic)
    {
        xi = x[ix];
        ai = a[ia];
        bi = b[ib];
        ci = c[ic];
        if_NA_dpq3_set(y[i], xi, ai, bi, ci)
        else
        {
            y[i] = f(xi, ai, bi, ci, i_1, i_2);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    FINISH_DPQ3;

    return sy;
}

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define DPQ3_1(A, FUN) dpq3_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN);
#define DPQ3_2(A, FUN) dpq3_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), FUN)

SEXP actuar_do_dpq3(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ3_1(args, dbetabinom);
    case  2:  return DPQ3_2(args, pbetabinom);
    case  3:  return DPQ3_2(args, qbetabinom);
    case  4:  return DPQ3_1(args, mbetabinom);
    case  5:  return DPQ3_1(args, dbetanbinom);
    case  6:  return DPQ3_2(args, pbetanbinom);
    case  7:  return DPQ3_2(args, qbetanbinom);
    case  8:  return DPQ3_1(args, mbetanbinom);
    default:
        error(_("internal error in actuar_do_dpq3"));
    }

    return args;                /* never used; to keep -Wall happy */
}
