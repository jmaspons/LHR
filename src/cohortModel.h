/*
 * Cohort model
 *
 */

int survdist(int n0, double surv, double varSurv, double limit, double ** pSurv);
int fertdist(double * years, int maxYear, int broods, int B, double surv, double varSurv, double * season, double ** pFert);

SEXP survdistR(SEXP args);
SEXP cohortModel(SEXP args);
SEXP exploreCohortModel(SEXP args);
