// TODO: nbinomialCompound and BetaNBinomialCompound
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> // probability distributions
#include "probability.h" // *betabinomial()


SEXP binomialCompound(SEXP args){
    SEXP ssize, sprob, slog, Res;
    int * log, maxSize, i, j;
    double * size, * prob, * res;
  
    args = CDR(args);
    PROTECT(ssize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sprob = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slog = coerceVector(CAR(args), LGLSXP));
    
    size = REAL(ssize);
    prob = REAL(sprob);
    log = LOGICAL(slog);
    
    maxSize = LENGTH(ssize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxSize));
    res = REAL(Res);
    memset(res, 0, maxSize * sizeof(double)); // initialise res
    // Rprintf("\tmaxX:%i\tprob=%.2f,\tlog=%i\n", maxSize, * prob, * log);
    for  (i=0; i < maxSize; i++){
      for (j=0; j <= i; j++){
        res[j] += dbinom(j, i, * prob, * log) * size[i];
        // Rprintf("\ti:%i; j=%i; p=%.6f; cum=%.6f\n", i, j, dbinom(j, i, * prob, * log) * size[i], res[j]);
      }
    }
    
    UNPROTECT(4);
    return (Res);
}


SEXP BetaBinomialCompound(SEXP args){
    SEXP ssize, sshape1, sshape2, slog, Res;
    int * log, maxSize, i, j;
    double * size, * shape1, * shape2, * res;
  
    args = CDR(args);
    PROTECT(ssize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sshape1 = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sshape2 = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slog = coerceVector(CAR(args), LGLSXP));
    
    size = REAL(ssize);
    shape1 = REAL(sshape1);
    shape2 = REAL(sshape2);
    log = LOGICAL(slog);
    
    maxSize = LENGTH(ssize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxSize));
    res = REAL(Res);
    memset(res, 0, maxSize * sizeof(double)); // initialise res
    // Rprintf("\tmaxX:%i\tshape1=%.2f; shape2=%.2f\tlog=%i\n", maxSize, * shape1, * shape2, * log);
    for  (i=0; i < maxSize; i++){
      for (j=0; j <= i; j++){
        res[j] += dbetabinom(j, i, * shape1, * shape2, * log) * size[i];
        // Rprintf("\ti:%i; j=%i; p=%.6f; cum=%.6f\n", i, j, dbetabinom(j, i, * shape1, * shape2, * log) * size[i], res[j]);
      }
    }
    
    UNPROTECT(5);
    return (Res);
}


//smax<- max(x) + max(y)
SEXP distrisum(SEXP args){
    SEXP sx, spx, sy, spy, smax, Res;
    int * x, * y, * max, xy, i, j;
    double * px, * py, * res;
    
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
    Res = PROTECT(allocVector(REALSXP, * max + 1));
    res = REAL(Res);
    memset(res, 0, (* max + 1) * sizeof(double)); // initialise res
    // Rprintf("\tmaxX:%i\tlenX=%i,lenY=%i\n", * max, LENGTH(sx), LENGTH(sy));
    for  (i=0; i < LENGTH(sx); i++){
      for (j=0; j < LENGTH(sy); j++){
        xy = x[i] + y[j];
        res[xy] += px[i] * py[j];
        // Rprintf("\tx:%i; p=%.4f; pp=%.4f\txx=%i; px=%.2f; xy=%i; py=%.2f\n", xy, res[xy], px[i] * py[j], x[i], px[i], y[j], py[j]);
      }
    }

    UNPROTECT(6);
    return (Res);
}
