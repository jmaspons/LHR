/*
 * AUTHOR
 *  Copyright (C) 2012, Joan Maspons, joanmaspons@gmail.com
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * DESCRIPTION
 *   TODO
 *   Simulates the demography of a cohort of n0 individuals whith life history
 *   parameters xxx.
 *   To compute the beta binomial probability, call dbetabinom(x, size, shape1, shape2).
 *   This checks for argument validity, and calls dbetabinom_raw().
 */

#include <R.h>
#include <Rmath.h>
#include "probability.h"
#include "dpq.h"

//TODO Pointers to functions cplusplus to use different distributions. Better overloaded functions because they have different parameters

int survdist(int n0, double surv, double varSurv, double limit, double ** pSurv)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n0) || ISNAN(surv) || ISNAN(limit)) {
//       pSurv = NULL;
        return 0;
    }
#endif

    if (surv < 0 || surv > 1 || n0 < 0 || varSurv < 0 || limit < 0) {
REprintf("invalid parameters.");

        return 0;
    }

    double * ptrSurv, maxYear, p;
    int i;


    if (varSurv == 0) {
        if (limit < 1 && limit > 0) { // limit = quantile omited
            maxYear = qnbinom(limit, n0, surv, /*l._t.*/ FALSE, /*log_p*/ FALSE);
        } else { // limit = max x
            maxYear = R_D_forceint(limit);
        }

        ptrSurv = Calloc(maxYear+1, double); //Allows to retun a vector

        if (ptrSurv == NULL) {
            REprintf("pSurv == NULL\n");
            * pSurv = ptrSurv;
            return 0;
        }

        for (i = 0; i <= maxYear; i++) {
            ptrSurv[i] = dnbinom(i, n0, surv, /*log_p*/ FALSE);
            p += ptrSurv[i];
        }

    } else { // varSurv
        double shape[2];
        fbeta(surv, varSurv, shape);
        if (ISNAN(shape[0]) || ISNAN(shape[1]) || shape[0] == 0 || shape[1] == 0) {
// Rprintf("surv=%.2f & var=%.2f out of the beta distribution domain", surv, varSurv);
            * pSurv = NULL;
            return 0;
        }

        if (limit < 1 && limit > 0) { // limit = quantile omited
            double minShape1 = 1, maxX = 15000; // it gave huge nsurv with limit = 0.05

            if (shape[0] < minShape1 && pbetanbinom(maxX, n0, shape[0], shape[1], /*l._t.*/ FALSE, /*log_p*/ FALSE) > limit){
                * pSurv = NULL;
                return -1;
            }
            maxYear = qbetanbinom(limit, n0, shape[0], shape[1], /*l._t.*/ FALSE, /*log_p*/ FALSE);
        } else { // limit = max x
            maxYear = R_D_forceint(limit);
        }

        ptrSurv = Calloc(maxYear + 1, double); //Allows to retun a vector
        if (ptrSurv == NULL) {
            REprintf("pSurv == NULL\n");
            * pSurv = ptrSurv;
            return 0;
        }

        for (i = 0; i <= maxYear; i++) {
            ptrSurv[i] = R_pow(M_E, dbetanbinom(i, n0, shape[0], shape[1], /*log_p*/ TRUE));
            p += ptrSurv[i];
        }
    }

    for (i = 0; i <= maxYear; i++) {
        ptrSurv[i] /= p; //correct for omited probability
    }

    * pSurv = ptrSurv;

    return maxYear;
}

int fertdist(double * years, int maxYear, int broods, int B, double surv, double varSurv, double * season, double ** pFert)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(* years) || ISNAN(maxYear)|| ISNAN(broods) || ISNAN(B) || ISNAN(surv) || ISNAN(varSurv)) {
//         pFert = NULL;
        return 0;
    }
#endif

    if (R_D_negInonint(maxYear) || R_D_negInonint(broods) || R_D_negInonint(B) || surv < 0 || surv > 1 || varSurv < 0) {
        REprintf("invalid parameters.");
//         pFert = NULL;
        return 0;
    }

    double * ptrFert;
    int i, j, maxFit;

    maxFit = maxYear * B;
    ptrFert = Calloc(maxFit + 1, double);
    if (ptrFert == NULL) {
        REprintf("ptrFert == NULL\n");
        * pFert = ptrFert;

        return 0;
    }

    if (season == NULL) {
        if (varSurv == 0) {
            for (i = 0; i <= maxFit; i++) {
                for (j = maxYear; i <= j * B; j--) { // i <= j * B
                    ptrFert[i] += years[j] * dbinom(i, j * B, surv, FALSE); //TODO Call to dbetabinom_raw
                }
            }
        } else { // varSurv
            double shape[2];
            fbeta(surv, varSurv, shape);
            if (ISNAN(shape[0]) || ISNAN(shape[1]) || shape[0] == 0 || shape[1] == 0) { //TODO:segfault
// Rprintf("surv=%.2f & var=%.2f out of the beta distribution domain\n", surv, varSurv);
                * pFert = ptrFert;

                return 0;
            }

            for (i = 0; i <= maxFit; i++) {
                for (j = maxYear; i <= j * B; j--) { // i <= j * B
                    ptrFert[i] += years[j] * R_pow(M_E, dbetabinom(i, j * B, shape[0], shape[1], TRUE)); //TODO Call to dbetabinom_raw
                }
            }
        }
    } else { // seasonality
        double tmp;
        int k, clutch[broods];
        tmp = B % broods;
        B -= tmp;

        for (i = 0; i <= broods; i++) {
            clutch[i] = (int) (B / broods);
            if (tmp > 0) {
                clutch[i] += 1;
                tmp -= 1;
            }
        }
        i = 0;
        while (tmp > 0) {
            clutch[i] += 1;
            tmp -= 1;
            i++;
            if (i > broods - 1) i = 0;
        }

        if (varSurv == 0) {
            for (i = 0; i <= maxFit; i++) {
                for (j = maxYear; i <= j * clutch[0]; j--) { // i <= j * clutch[k]       clutch[0] == max(clutch)
                    for (k = 0; k < broods; k++) {
                        ptrFert[i] += years[j] * dbinom(i, j * clutch[k], surv * season[k], FALSE) / broods; //TODO Call to dbetabinom_raw
                    }
                }
            }

        } else { // seasonality + varSurv
            double shape[broods][2];
            for (i = 0; i < broods; i++) {
                fbeta(surv * season[i], varSurv, shape[i]);

                if (ISNAN(shape[i][0]) || ISNAN(shape[i][1]) || shape[i][0] == 0 || shape[i][1] == 0) { //TODO: R_NaN not recognized? TODO:segfault
//Rprintf("surv * season=%.2f & var=%.2f out of the beta distribution domain\n", surv * season[i], varSurv);
                    * pFert = ptrFert;

                    return 0;
                }
            }

            for (i = 0; i <= maxFit; i++) {
                for (j = maxYear; i <= j * clutch[0]; j--) { // i <= j * clutch[k]       clutch[0] == max(clutch)
                    for (k = 0; k < broods; k++) {
                        if (i <= j * clutch[k]) // FIXME exp(R_NegInf) == 1 instead of 0
                            ptrFert[i] += years[j] * R_pow(M_E, dbetabinom(i, j * clutch[k], shape[k][0], shape[k][1], TRUE)) / broods; //TODO Call to dbetabinom_raw
                    }
                }
            }
        }
    }

    * pFert = ptrFert;

    return maxFit;
}


