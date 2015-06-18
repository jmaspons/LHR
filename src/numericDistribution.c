// TODO: nbinomialCompound and BetaNBinomialCompound
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> // probability distributions
#include "probability.h" // *betabinomial()


SEXP binomialCompound(SEXP args){
    SEXP ssize, spsize, sprob, slog, smaxsize, Res;
    int log, maxsize, i, j;
    double * size, * psize, prob, * res;
  
    args = CDR(args);
    PROTECT(ssize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(spsize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sprob = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slog = coerceVector(CAR(args), LGLSXP));
    args = CDR(args);
    PROTECT(smaxsize = coerceVector(CAR(args), REALSXP));
    
    size = REAL(ssize);
    psize = REAL(spsize);
    prob = asReal(sprob);
    log = asLogical(slog);
    maxsize = asInteger(smaxsize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxsize + 1));
    res = REAL(Res);
    memset(res, 0, (maxsize + 1) * sizeof(double)); // initialise res
    // Rprintf("\tmaxX:%i\tprob=%.2f,\tlog=%i\n", maxsize, prob, log);
    for  (i=0; i < LENGTH(ssize); i++){
      for (j=0; j <= size[i]; j++){
          res[j] += dbinom((double) j, size[i], prob, log) * psize[i];
        // Rprintf("\ti:%i; x=%i; p*psize=%.6f; xcum=%.6f; p=%.3f, size=%.0f, psize=%.2f\n", i, j, dbinom((double) j, size[i], prob, log) * psize[i], res[j], dbinom((double) j, size[i], prob, log), size[i], psize[i]);
      }
    }
    
    UNPROTECT(6);
    return (Res);
}


SEXP BetaBinomialCompound(SEXP args){
    SEXP ssize, spsize, sshape1, sshape2, slog, smaxsize, Res;
    int log, maxsize, i, j;
    double * size, * psize, * shape1, * shape2, * res;
  
    args = CDR(args);
    PROTECT(ssize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(spsize = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sshape1 = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(sshape2 = coerceVector(CAR(args), REALSXP));
    args = CDR(args);
    PROTECT(slog = coerceVector(CAR(args), LGLSXP));
    args = CDR(args);
    PROTECT(smaxsize = coerceVector(CAR(args), REALSXP));
    
    size = REAL(ssize);
    psize = REAL(spsize);
    shape1 = REAL(sshape1);
    shape2 = REAL(sshape2);
    log = asLogical(slog);
    maxsize = asInteger(smaxsize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxsize + 1));
    res = REAL(Res);
    memset(res, 0, (maxsize + 1) * sizeof(double)); // initialise res
    // Rprintf("\tmaxX:%i\tshape1=%.2f; shape2=%.2f\tlog=%i\n", maxSize, * shape1, * shape2, * log);
    for  (i=0; i < LENGTH(ssize); i++){
      for (j=0; j <= size[i]; j++){
        res[j] += dbetabinom(j, size[i], * shape1, * shape2, log) * psize[i];
        // Rprintf("\ti:%i; j=%i; p=%.6f; cum=%.6f\n", i, j, dbetabinom(j, size[i], * shape1, * shape2, log) * psize[i], res[j]);
      }
    }
    
    UNPROTECT(7);
    return (Res);
}


//smax<- max(x) + max(y)
SEXP distrisum(SEXP args){
    SEXP sx, spx, sy, spy, smax, Res;
    int * x, * y, max, xy, i, j;
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
    max = asInteger(smax);
    
    // res contain the probabilities of resulting distribution for values from 0 to smax
    Res = PROTECT(allocVector(REALSXP, max + 1));
    res = REAL(Res);
    memset(res, 0, (max + 1) * sizeof(double)); // initialise res
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
