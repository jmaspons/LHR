// TODO: nbinomialCompound and BetaNBinomialCompound
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> // probability distributions
#include "probability.h" // *betabinomial()


SEXP binomialCompound(SEXP args){
    SEXP ssize, spsize, sprob, slog, smaxsize, Res;
    int logP, maxsize, i, j;
    double * size, * psize, prob, p, * res;
  
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
    logP = asLogical(slog);
    maxsize = asInteger(smaxsize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxsize + 1));
    res = REAL(Res);
    memset(res, 0, (maxsize + 1) * sizeof(double)); // initialise res. Take care later about the logP scale
    // Rprintf("\tmaxX:%i\tprob=%.2f,\tlog=%i\n", maxsize, prob, logP);
    if (logP){
      for  (i=0; i < LENGTH(ssize); i++){
        for (j=0; j <= size[i]; j++){
          p = dbinom((double) j, size[i], prob, logP);
          // Rprintf("\ti:%i; x=%i; size=%.0f; psize=%.2f; p=%.3f; ## ", i, j, size[i], psize[i], p);
          // Rprintf("p + psize=%.6f; pcum=%.6f; pcumC=%6f; logpcumC=%.3f\n", p + psize[i], logspace_add(res[j], p + psize[i]), exp(res[j]) + exp(p+psize[i]), log(exp(res[j]) + exp(p+psize[i])));
          res[j] = logspace_add(res[j], p + psize[i]);
        }
      }
      // substract 0 for the initialization on the log scale
      for (i=0; i <= maxsize; i++){
        res[i] = logspace_sub(res[i], 0);
      }
    }else{
      for  (i=0; i < LENGTH(ssize); i++){
        for (j=0; j <= size[i]; j++){
          p = dbinom((double) j, size[i], prob, logP);
          res[j] += p * psize[i];
          // Rprintf("\ti:%i; x=%i; size=%.0f; psize=%.2f; p=%.3f; ## ", i, j, size[i], psize[i], p);
          // Rprintf("p * psize=%.6f; xcum=%.6f\n", p * psize[i], res[j]);
        }
      }
    }

    
    UNPROTECT(6);
    return (Res);
}


SEXP BetaBinomialCompound(SEXP args){
    SEXP ssize, spsize, sshape1, sshape2, slog, smaxsize, Res;
    int logP, maxsize, i, j;
    double * size, * psize, * shape1, * shape2, p, * res;
  
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
    logP = asLogical(slog);
    maxsize = asInteger(smaxsize);
    
    // res contain the probabilities of resulting distribution for values from 0 to maxSize
    Res = PROTECT(allocVector(REALSXP, maxsize + 1));
    res = REAL(Res);
    memset(res, 0, (maxsize + 1) * sizeof(double)); // initialise res. Take care later about the logP scale
    // Rprintf("\tmaxX:%i\tshape1=%.2f; shape2=%.2f\tlog=%i\n", maxsize, * shape1, * shape2, * logP);
    if (logP){
      for  (i=0; i < LENGTH(ssize); i++){
        for (j=0; j <= size[i]; j++){
          p = dbetabinom(j, size[i], * shape1, * shape2, logP);
          res[j] = logspace_add(res[j], p + psize[i]);
          // Rprintf("\ti:%i; x=%i; size=%.0f; psize=%.2f; p=%.3f; ## ", i, j, size[i], psize[i], p);
          // Rprintf("p + psize=%.6f; pcum=%.6f; pcumC=%6f; logpcumC=%.3f\n", p + psize[i], logspace_add(res[j], p + psize[i]), exp(res[j]) + exp(p+psize[i]), log(exp(res[j]) + exp(p+psize[i])));
        }
      }
      // substract 0 for the initialization on the log scale
      for (i=0; i <= maxsize; i++){
        res[i] = logspace_sub(res[i], 0);
      }
    }else{
      for  (i=0; i < LENGTH(ssize); i++){
        for (j=0; j <= size[i]; j++){
          p = dbetabinom(j, size[i], * shape1, * shape2, logP);
          res[j] += p * psize[i];
          // Rprintf("\ti:%i; x=%i; size=%.0f; psize=%.2f; p=%.3f; ## ", i, j, size[i], psize[i], p);
          // Rprintf("p + psize=%.6f; pcum=%.6f; pcumC=%6f; logpcumC=%.3f\n", p + psize[i], logspace_add(res[j], p + psize[i]), exp(res[j]) + exp(p+psize[i]), log(exp(res[j]) + exp(p+psize[i])));
        }
      }
    }
    
    UNPROTECT(7);
    return (Res);
}


//smax<- max(x) + max(y)
SEXP distrisum(SEXP args){
    SEXP sx, spx, sy, spy, slog, smax, Res;
    int * x, * y, max, xy, logP, i, j;
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
    PROTECT(slog = coerceVector(CAR(args), LGLSXP));
    args = CDR(args);
    PROTECT(smax = coerceVector(CAR(args), INTSXP));
    
    x = INTEGER(sx);
    px = REAL(spx);
    y = INTEGER(sy);
    py = REAL(spy);
    logP = asLogical(slog);
    max = asInteger(smax);
    
    // res contain the probabilities of resulting distribution for values from 0 to smax
    Res = PROTECT(allocVector(REALSXP, max + 1));
    res = REAL(Res);
    memset(res, 0, (max + 1) * sizeof(double)); // initialise res. Take care later about the logP scale
    // Rprintf("\tmaxX:%i\tlenX=%i,lenY=%i\n", * max, LENGTH(sx), LENGTH(sy));
    if (logP){
      for  (i=0; i < LENGTH(sx); i++){
        for (j=0; j < LENGTH(sy); j++){
          xy = x[i] + y[j];
          res[xy] = logspace_add(res[xy], px[i] + py[j]);
          // Rprintf("\tx:%i; p=%.4f; pp=%.4f\txx=%i; px=%.2f; xy=%i; py=%.2f\n", xy, res[xy], px[i] * py[j], x[i], px[i], y[j], py[j]);
        }
      }
      // substract 0 for the initialization on the log scale
      for (i=0; i <= max; i++){
        res[i] = logspace_sub(res[i], 0);
      }
    }else{
      for  (i=0; i < LENGTH(sx); i++){
        for (j=0; j < LENGTH(sy); j++){
          xy = x[i] + y[j];
          res[xy] += px[i] * py[j];
          // Rprintf("\tx:%i; p=%.4f; pp=%.4f\txx=%i; px=%.2f; xy=%i; py=%.2f\n", xy, res[xy], px[i] * py[j], x[i], px[i], y[j], py[j]);
        }
      }
    }
    
    UNPROTECT(7);
    return (Res);
}
