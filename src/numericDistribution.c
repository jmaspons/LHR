#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

//smax<- max(x) + max(y)
SEXP distrisum(SEXP args){
    SEXP sx, spx, sy, spy, smax, res;
    int * x, * y, * max, xy, i, j;
    double * px, * py, * pxy;
  
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
    PROTECT(res = allocVector(REALSXP, * max + 1));
    // Rprintf("\tmaxX:%i\tlenX=%i,lenY=%i\n", * max, LENGTH(sx), LENGTH(sy));
    for  (i=0; i < LENGTH(sx); i++){
      for (j=0; j < LENGTH(sy); j++){
        xy = x[i] + y[j];
        REAL(res)[xy] += px[i] * py[j];
       // Rprintf("\tx:%i; xx=%i; xy=%i; \tp=%.2f\n", xy, x[i], y[j], px[i] * py[j]);
      }
    }
    UNPROTECT(6);
    return (res);
}