#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "cohortModel.h"
#include "environment.h"


SEXP cohortModel(SEXP args)
{
    SEXP sn0, ssurvA, svarSurvA, slimit, sbroods, sB, ssurvJ, svarSurvJ, sseason, sfitdist;
    int * n0, * B, * broods, n, i;
    double * survA, * survJ, * varSurvA, * varSurvJ, * limit, * season, * pSurv, * fitdist;

    args = CDR(args);
    PROTECT(sn0 = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(ssurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(svarSurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slimit = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sbroods = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(sB = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(ssurvJ = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(svarSurvJ = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sseason = coerceVector(CAR(args), REALSXP));

    n0 = INTEGER(sn0);
    survA = REAL(ssurvA);
    varSurvA = REAL(svarSurvA);
    limit = REAL(slimit);
    broods = INTEGER(sbroods);
    B = INTEGER(sB);
    survJ = REAL(ssurvJ);
    varSurvJ = REAL(svarSurvJ);
    season = REAL(sseason);

// Rprintf("n0=%i, survA=%.2f, varSurvA=%.2f, limit=%.2f, broods=%i, B=%i, survJ=%.2f, varSurvJ=%.2f, season[0]=%.2f\n", * n0, * survA, * varSurvA, * limit, * broods, * B, * survJ, * varSurvJ, * season);

    n = survdist(* n0, * survA, * varSurvA, * limit, & pSurv);
// Rprintf("nSurvdist:%i", n);
    if (n > 0) {
        n = fertdist(pSurv, n, * broods, * B, * survJ, * varSurvJ, season, & fitdist);
// Rprintf("\tnFertdist:%i\n", n);
    }

    if (n > 0 ) {
        PROTECT(sfitdist = allocVector(REALSXP, n + 1));

        for (i = 0; i <= n; i++) {
            REAL(sfitdist)[i] = fitdist[i];
        }

    } else {
        PROTECT(sfitdist = allocVector(REALSXP, 1));
        REAL(sfitdist)[0] = R_NaN;
    }

    Free(pSurv);
    Free(fitdist);
    UNPROTECT(10);

    return sfitdist;
}

SEXP survdistR(SEXP args)
{
    SEXP sn0, ssurvA, svarSurvA, slimit, spSurv;
    int * n0, n, i;
    double * survA, * varSurvA, * limit, *pSurv;

    args = CDR(args);
    PROTECT(sn0 = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(ssurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(svarSurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slimit = coerceVector(CAR(args), REALSXP));

    n0 = INTEGER(sn0);
    survA = REAL(ssurvA);
    varSurvA = REAL(svarSurvA);
    limit = REAL(slimit);

    // Rprintf("n0=%i, survA=%.2f, varSurvA=%.2f, limit=%.2f\n", * n0, * survA, * varSurvA, * limit);

    n = survdist(* n0, * survA, * varSurvA, * limit, & pSurv);
    if (n > 0 ) {
        PROTECT(spSurv = allocVector(REALSXP, n + 1));

        for (i = 0; i <= n; i++) {
            REAL(spSurv)[i] = pSurv[i];
        }

    } else {
        PROTECT(spSurv = allocVector(REALSXP, 1));
        REAL(spSurv)[0] = R_NaN;
    }

    Free(pSurv);
    UNPROTECT(5);

    return spSurv;
}

SEXP exploreCohortModel(SEXP args)
{
    SEXP sn0, ssurvA, svarSurvA, slimit, sbroods, sB, ssurvJ, svarSurvJ, sseasonAmpl, sbreedInterval, res, dimnames, colnames;
    int * n0, * B, * broods, i, j, nsurv, nfit, ncol, nrow, found, interannual;
    double * survA, * survJ, * varSurvA, * varSurvJ, * limit, * seasonAmpl, * breedInterval, * seasonPattern, * pSurv, * fitdist, * resPtr, cumsum, replace;

    args = CDR(args);
    PROTECT(sn0 = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(ssurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(svarSurvA = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slimit = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sbroods = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(sB = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(ssurvJ = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(svarSurvJ = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sseasonAmpl = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sbreedInterval = coerceVector(CAR(args), REALSXP));

    nrow = LENGTH(sn0);
    ncol = 5; // mean, var, G, pReplacement, error
    replace = 2;

    PROTECT(res = allocMatrix(REALSXP, nrow, ncol)); // matrix increase index of the first dimension faster

    n0 = INTEGER(sn0);
    survA = REAL(ssurvA);
    varSurvA = REAL(svarSurvA);
    limit = REAL(slimit);
    broods = INTEGER(sbroods);
    B = INTEGER(sB);
    survJ = REAL(ssurvJ);
    varSurvJ = REAL(svarSurvJ);
    seasonAmpl = REAL(sseasonAmpl);
    breedInterval = REAL(sbreedInterval);

    resPtr = REAL(res);

    Rprintf("case 1 / %i\tn0=%i, sA=%.2f, varSa=%.2f  ### B=%i, sJ=%.2F, varSj=%.2f, ampl=%.2f\n", nrow, n0[0], survA[0], varSurvA[0], B[0], survJ[0], varSurvJ[0], seasonAmpl[0]);
    nsurv = survdist(* n0, * survA, * varSurvA, * limit, & pSurv);
// Rprintf("Survdist done. Nsurv=%i\n", nsurv);
    if (nsurv > 0) {
        interannual = 0;

        if (* seasonAmpl > 0) {
            if (seasonalPattern( * broods, * breedInterval, * seasonAmpl, & seasonPattern) < 0)
                interannual = 1;
        } else {
            seasonPattern = NULL;
        }

        nfit = fertdist(pSurv, nsurv, * broods, * B, * survJ, * varSurvJ, seasonPattern, & fitdist);
// Rprintf("Fertdist done. Nfit=%i\n", nfit);
        if (nfit > 0) {
            if (interannual == 1) {
                resPtr[4 * nrow] = - * broods;
            } else {
                resPtr[4 * nrow] = 0;
            }
            for (j = 0, found = 0, cumsum = 0; j <= nfit; j++) {
                resPtr[0] += fitdist[j] * j / * n0; // mean = sum(P(i) * i / n0))
                if (!found) {
                    cumsum += fitdist[j];
                    if (j / * n0 >= replace) {
                        resPtr[3 * nrow] = 1 - cumsum; // pReplace = P(i) | i == 2 * n0
                        found = 1;
                    }
                }
            }

            for (j = 0; j <= nfit; j++) {
                resPtr[1 * nrow] += fitdist[j] * R_pow_di(j / * n0 - resPtr[0], 2); // var = sum(P(i) * (i / n0 - mean)^2)
            }

            resPtr[2 * nrow] = resPtr[0] - 2 * resPtr[1 * nrow] / resPtr[0]; // G = mean - 2 * var / mean
        }

    } else {
        nfit = 0;
        seasonPattern = NULL;
        pSurv = NULL;
        fitdist = NULL;
    }

    if (nfit == 0) {
        fitdist = NULL;
        for (j = 0; j < 4; j++) {
            resPtr[j * nrow] = R_NaReal;
        }
        resPtr[4 * nrow] = nsurv > 0 ? 2 : (nsurv == 0? 1 : 1.5); // ERROR 1: survdist, 2: fitdist, +.5: maxYear >15000, -broods: interannual breeding pattern
    }

    for (i = 1; i < nrow; i++) {
        R_CheckUserInterrupt();

        Rprintf("case %i / %i\tn0=%i, sA=%.2f, varSa=%.2f  ### B=%i, sJ=%.2F, varSj=%.2f, ampl=%.2f\n", i+1, nrow,  n0[i], survA[i], varSurvA[i], B[i], survJ[i], varSurvJ[i], seasonAmpl[i]);

        if (n0[i] != n0[i-1] || survA[i] != survA[i-1] || varSurvA[i] != varSurvA[i-1] || limit[i] != limit[i-1]) {
            Free(pSurv);
            nsurv = survdist(n0[i], survA[i], varSurvA[i], limit[i], & pSurv);
        }
        Rprintf("\tnsurv=%i", nsurv);
        if (nsurv > 0) {
            Free(fitdist);
            interannual = 0;

            if (seasonAmpl[i] > 0) {
                if (broods[i] != broods[i-1] || breedInterval[i-1] != breedInterval[i-1] || seasonAmpl[i] != seasonAmpl[i-1]) {
                    Free(seasonPattern);
                    if (seasonalPattern(broods[i], breedInterval[i], seasonAmpl[i], & seasonPattern) < 0)
                        interannual = 1;
                }
            } else {
                seasonPattern = NULL;
            }

            nfit = fertdist(pSurv, nsurv, broods[i], B[i], survJ[i], varSurvJ[i], seasonPattern, & fitdist);
            Rprintf(", nfit=%i", nfit);
            if (nfit > 0) {
                if (interannual == 1) {
                    resPtr[4 * nrow + i] = -broods[i];
                } else {
                    resPtr[4 * nrow + i] = 0;
                }

                for (j = 0, found = 0, cumsum = 0; j <= nfit; j++) {
                    resPtr[0 * nrow + i] += fitdist[j] * j / n0[i]; // mean = sum(P(i) * i / n0))
                    if (!found) {
                        cumsum += fitdist[j];
                        if (j / n0[i] >= replace) {
                            resPtr[3 * nrow + i] = 1 - cumsum; // pReplace = P(i) | i == 2 * n0
                            found = 1;
                        }
                    }
                }

                for (j = 0; j <= nfit; j++) {
                    resPtr[1 * nrow + i] += fitdist[j] * R_pow_di(j / n0[i] - resPtr[0 * nrow + i], 2); // var = sum(P(i) * (i / n0 - mean)^2)
                }

                resPtr[2 * nrow + i] = resPtr[0 * nrow + i] - 2 * resPtr[1 * nrow + i] / resPtr[0 * nrow + i]; // G = mean - 2 * var / mean
            }

        } else {
            nfit = 0;
            Rprintf(", nfit=%i", nfit);
        }

        if (nfit == 0) {
            for (j = 0; j < 4; j++) {
                resPtr[j * nrow + i] = R_NaReal;
            }
            resPtr[4 * nrow + i] = nsurv > 0 ? 2 : (nsurv == 0? 1 : 1.5); // ERROR 1: survdist, 2: fitdist, +.5: maxYear >15000, -broods: interannual breeding pattern
        }
        Rprintf("\tmean=%.2f, var=%.2f, G=%.2f, Prep=%.2f, error=%.0f\n", resPtr[0 * nrow + i], resPtr[1 * nrow + i], resPtr[2 * nrow + i], resPtr[3 * nrow + i], resPtr[4 * nrow + i]);
    }
    Rprintf("Simulation done!\n");
    Free(seasonPattern);
    Free(pSurv);
    Free(fitdist);
    //Pàg. 110 i 125 R_exts.pdf
//     PROTECT(colnames = allocVector(VECSXP, 5));
    //     SET_VECTOR_ELT(colnames, 0, "mean");
    //     SET_VECTOR_ELT(colnames, 1, "var");
    //     SET_VECTOR_ELT(colnames, 2, "G");
    //     SET_VECTOR_ELT(colnames, 3, "Preplace");
    //     SET_VECTOR_ELT(colnames, 4, "error");

//     PROTECT(dimnames = allocVector(VECSXP, 2));
//     SET_VECTOR_ELT(dimnames, 1, colnames);
//     dimnamesgets(gradient, dimnames);
//     setAttrib(res, R_DimNamesSymbol, dimnames);

//     PROTECT(dimnames = allocVector(VECSXP, 2));
//     SET_VECTOR_ELT(dimnames, 0, getAttrib(x, R_NamesSymbol));
//     SET_VECTOR_ELT(dimnames, 1, getAttrib(y, R_NamesSymbol));
//     setAttrib(ans, R_DimNamesSymbol, dimnames);
//     for (i=0; i<nrow; i++){
// Rprintf("[%i] mu=%.2f, var=%.2f, G=%.2f, pRep=%.2f, err=%.0f\n", i, REAL(res)[0 * nrow + i] , REAL(res)[1 * nrow + i], REAL(res)[2 * nrow + i], REAL(res)[3 * nrow + i], REAL(res)[4 * nrow + i]);
//     }

    UNPROTECT(11);

    return (res);
}
// #include <R_ext/PrtUtil.h>
// SEXP showArgs(SEXP args)
// {
//   int i;
//   args = CDR(args); /* skip ’name’ */
//   for(i = 0; args != R_NilValue; i++, args = CDR(args)) {
//     const char *name =
//     isNull(TAG(args)) ? "" : CHAR(PRINTNAME(TAG(args)));
//     SEXP el = CAR(args);
//     if (length(el) != 0) {
//       Rprintf("[%d] ’%s’ R type, length 0\n", i+1, name);
//     continue;
//     }
//     switch(TYPEOF(el)) {
//       case REALSXP:
//         Rprintf("[%d] ’%s’ %f\n", i+1, name, REAL(el)[0]);
//         break;
//       case LGLSXP:
//       case INTSXP:
//         Rprintf("[%d] ’%s’ %d\n", i+1, name, INTEGER(el)[0]);
//         break;
//       case CPLXSXP:
//       {
//         Rcomplex cpl = COMPLEX(el)[0];
//         Rprintf("[%d] ’%s’ %f + %fi\n", i+1, name, cpl.r, cpl.i);
//       }
//         break;
//       case STRSXP:
//         Rprintf("[%d] ’%s’ %s\n", i+1, name,
//         CHAR(STRING_ELT(el, 0)));
//         break;
//       default:
//         Rprintf("[%d] ’%s’ R type\n", i+1, name);
//     }
//   }
//   return(R_NilValue);
// }

/*
// R-exts.pdf pg. 117 .External
#include <R.h>
#include <Rinternals.h>
SEXP convolveE(SEXP args)
{
  R_len_t i, j, na, nb, nab;
  double *xa, *xb, *xab;
  SEXP a, b, ab;
  PROTECT(a = coerceVector(CADR(args), REALSXP));
  PROTECT(b = coerceVector(CADDR(args), REALSXP));
  ...
}
first = CADR(args);
second = CADDR(args);
third = CADDDR(args);
fourth = CAD4R(args);
//provide convenient ways to access the first four arguments. More generally we can use the CDR and CAR macros as in:
args = CDR(args); a = CAR(args);
args = CDR(args); b = CAR(args);


// R-exts.pdf pg. 116 .Call
#include <R.h>
#include <Rdefines.h>
SEXP convolve2(SEXP a, SEXP b)
{
  R_len_t i, j, na, nb, nab;
  double *xa, *xb, *xab;
  SEXP ab;

PROTECT(a = AS_NUMERIC(a));
PROTECT(b = AS_NUMERIC(b));
na = LENGTH(a); nb = LENGTH(b); nab = na + nb - 1;
PROTECT(ab = NEW_NUMERIC(nab));
xa = NUMERIC_POINTER(a); xb = NUMERIC_POINTER(b);
xab = NUMERIC_POINTER(ab);
for(i = 0; i < nab; i++) xab[i] = 0.0;
for(i = 0; i < na; i++)
  for(j = 0; j < nb; j++) xab[i + j] += xa[i] * xb[j];
  UNPROTECT(3);
return(ab);
}

#include <R.h>
#include <Rinternals.h>
SEXP convolve2(SEXP a, SEXP b)
{
  R_len_t i, j, na, nb, nab;
  double *xa, *xb, *xab;
  SEXP ab;
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  na = length(a); nb = length(b); nab = na + nb - 1;
  PROTECT(ab = allocVector(REALSXP, nab));
  xa = REAL(a); xb = REAL(b);
  xab = REAL(ab);
  for(i = 0; i < nab; i++) xab[i] = 0.0;
  for(i = 0; i < na; i++)
    for(j = 0; j < nb; j++) xab[i + j] += xa[i] * xb[j];
    UNPROTECT(3);
  return(ab);
}

###########

// R-exts.pdf pg. 107 Handling R objects
SEXP ab;
....
PROTECT(ab = allocVector(REALSXP, 2));
REAL(ab)[0] = 123.45;
REAL(ab)[1] = 67.89;
UNPROTECT(1);

//and then those in ‘Rdefines.h’:

SEXP ab;
....
PROTECT(ab = NEW_NUMERIC(2));
NUMERIC_POINTER(ab)[0] = 123.45;
NUMERIC_POINTER(ab)[1] = 67.89;
UNPROTECT(1);
*/
