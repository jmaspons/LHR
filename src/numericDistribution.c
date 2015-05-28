#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

//smax<- max(x) + max(y)
SEXP distrisum(SEXP args){
    SEXP sx, spx, sy, spy, smax, Res;
    int * x, * y, * max, xy, i, j;
    double * px, * py, * res; //, * pxy;
  
    args = CDR(args);
    PROTECT(sx = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(spx = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sy = coerceVector(CAR(args), INTSXP));
    args = CDR(args);
    PROTECT(spy = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(smax = coerceVector(CAR(args), INTSXP));
    
    x = INTEGER(sx);
    px = REAL(spx);
    y = INTEGER(sy);
    py = REAL(spy);
    max = INTEGER(smax);
    
    // res contain the probabilities of resulting distribution for values from 0 to smax
    res = Calloc(* max + 1, double);
    // Rprintf("\tmaxX:%i\tlenX=%i,lenY=%i\n", * max, LENGTH(sx), LENGTH(sy));
    for  (i=0; i < LENGTH(sx); i++){
      for (j=0; j < LENGTH(sy); j++){
        xy = x[i] + y[j];
        res[xy] += px[i] * py[j];
        // Rprintf("\tx:%i; p=%.4f; pp=%.4f\txx=%i; px=%.2f; xy=%i; py=%.2f\n", xy, res[xy], px[i] * py[j], x[i], px[i], y[j], py[j]);
      }
    }
    
    PROTECT(Res = allocVector(REALSXP, * max + 1));
    for (i=0; i <= * max; i++){
      REAL(Res)[i] = res[i];
    }
    
    Free(res);
    UNPROTECT(6);
    return (Res);
}