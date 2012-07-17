/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000 and Feb, 2001.
 *
 *    dbetanbinom_mu(): Martin Maechler, June 2008
 *
 *  Merge in to R:
 *	Copyright (C) 2000--2008, The R Core Team
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
 *
 *   Computes the negative binomial distribution. For integer n,
 *   this is probability of x failures before the nth success in a
 *   sequence of Bernoulli trials. We do not enforce integer n, since
 *   the distribution is well defined for non-integers,
 *   and this can be useful for e.g. overdispersed discrete survival times.
 */

#include <R.h>
#include <Rmath.h>
#include "betaBinom.h"
#include "dpq.h"
#include "R_ext/Visibility.h"
//TODO
double attribute_hidden dbetanbinom_raw(double x, double size, double shape1, double shape2, int give_log)
{
    double ans;
 
    if (x < 0 || !R_FINITE(x)) return R_D__0;
    
    ans = pochhammer(shape1, size) * pochhammer(size, x) * pochhammer(shape2, x) / (gamma(x+1) * pochhammer(shape1 + shape2, size) * pochhammer(size + shape1 + shape2, x)); //factorial(x) = gamma(x+1)

    return R_D_qIv(ans);
}

double dbetanbinom(double x, double size, double shape1, double shape2, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
        return x + size + shape1 + shape2;
#endif

    if (shape1 <= 0 || shape2 <= 0 || size < 0) return R_NaN; //TODO  R_D_negInonint(size)
    R_D_nonint_check(x);
    x = R_D_forceint(x);

    return dbetanbinom_raw(x, size, shape1, shape2, give_log);
}

// double dnbinom_mu(double x, double size, double mu, int give_log)
// {
//     /* originally, just set  prob :=  size / (size + mu)  and called dbinom_raw(),
//      * but that suffers from cancellation when   mu << size  */
//     double ans, p;
// 
// #ifdef IEEE_754
//     if (ISNAN(x) || ISNAN(size) || ISNAN(mu))
//         return x + size + mu;
// #endif
// 
//     if (mu < 0 || size < 0) ML_ERR_return_NAN;
//     R_D_nonint_check(x);
//     if (x < 0 || !R_FINITE(x)) return R_D__0;
//     x = R_D_forceint(x);
//     if(x == 0)/* be accurate, both for n << mu, and n >> mu :*/
// 	return R_D_exp(size * (size < mu ? log(size/(size+mu)) : log1p(- mu/(size+mu))));
//     if(x < 1e-10 * size) { /* don't use dbinom_raw() but MM's formula: */
// 	/* FIXME --- 1e-8 shows problem; rather use algdiv() from ./toms708.c */
// 	return R_D_exp(x * log(size*mu / (size+mu)) - mu - lgamma(x+1) +
// 		       log1p(x*(x-1)/(2*size)));
//     }
//     /* else: no unnecessary cancellation inside dbinom_raw, when
//      * x_ = size and n_ = x+size are so close that n_ - x_ loses accuracy
//      */
//     ans = dbinom_raw(size, x+size, size/(size+mu), mu/(size+mu), give_log);
//     p = ((double)size)/(size+x);
//     return((give_log) ? log(p) + ans : p * ans);
// }

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-8 The R Core Team
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
 *  DESCRIPTION
 *
 *	The distribution function of the negative binomial distribution.
 *
 *  NOTES
 *
 *	x = the number of failures before the n-th success
 */

double attribute_hidden pbetanbinom_raw(double x, double size, double shape1, double shape2, int log_p)
{
    int i;
    double ans = 0;
    for (i = 0; i <= x; ++i) {
        ans += dbetanbinom_raw(i, size, shape1, shape2, log_p);
    }
    return ans;
}

double pbetanbinom(double x, double size, double shape1, double shape2, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
	return x + size + shape1 + shape2;
    if(!R_FINITE(size) || !R_FINITE(shape1) || !R_FINITE(shape2)) return R_NaN;
#endif
    if (size <= 0 || shape1 <= 0 || shape2 <= 0) return R_NaN;

    if (x < 0) return R_DT_0;
    if (!R_FINITE(x)) return R_DT_1;
    x = floor(x + 1e-7);
    
    double ans = pbetanbinom_raw(x, size, shape1, shape2, log_p);
    if (!lower_tail) ans = 1 - ans;
    return ans;
}

// double pbetanbinom_mu(double x, double size, double mu, int lower_tail, int log_p)
// {
// #ifdef IEEE_754
//     if (ISNAN(x) || ISNAN(size) || ISNAN(mu))
// 	return x + size + mu;
//     if(!R_FINITE(size) || !R_FINITE(mu))	ML_ERR_return_NAN;
// #endif
//     if (size <= 0 || mu < 0)	ML_ERR_return_NAN;
// 
//     if (x < 0) return R_DT_0;
//     if (!R_FINITE(x)) return R_DT_1;
//     x = floor(x + 1e-7);
//     /* return
//      * pbeta(pr, size, x + 1, lower_tail, log_p);  pr = size/(size + mu), 1-pr = mu/(size+mu)
//      *
//      *= pbeta_raw(pr, size, x + 1, lower_tail, log_p)
//      *            x.  pin   qin
//      *=  bratio (pin,  qin, x., 1-x., &w, &wc, &ierr, log_p),  and return w or wc ..
//      *=  bratio (size, x+1, pr, 1-pr, &w, &wc, &ierr, log_p) */
//     {
// 	int ierr;
// 	double w, wc;
// 	bratio(size, x+1, size/(size+mu), mu/(size+mu), &w, &wc, &ierr, log_p);
// 	if(ierr)
// 	    MATHLIB_WARNING(_("pbetanbinom_mu() -> bratio() gave error code %d"), ierr);
// 	return lower_tail ? w : wc;
//     }
// }

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2008 The R Core Team
 *  Copyright (C) 2005 The R Foundation
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
 *  SYNOPSIS
 *
 *	#include <Rmath.h>
 *	double qbetanbinom(double p, double size, double shape1, double shape2,
 *                     int lower_tail, int log_p)
 *
 *  DESCRIPTION
 *
 *	The quantile function of the negative binomial distribution.
 *
 *  NOTES
 *
 *	x = the number of failures before the n-th success
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 */

static double
do_search(double y, double *z, double p, double n, double shape1, double shape2, double incr)
{
    if(*z >= p) {
			/* search to the left */
	for(;;) {
	    if(y == 0 ||
	       (*z = pbetanbinom(y - incr, n, shape1, shape2, /*l._t.*/TRUE, /*log_p*/FALSE)) < p)
		return y;
	    y = fmax2(0, y - incr);
	}
    }
    else {		/* search to the right */

	for(;;) {
	    y = y + incr;
	    if((*z = pbetanbinom(y, n, shape1, shape2, /*l._t.*/TRUE, /*log_p*/FALSE)) >= p)
		return y;
	}
    }
}


double qbetanbinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p)
{
    double P, Q, mu, sigma, gamma, z, y;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
	return p + size + shape1 + shape2;
#endif
    if (shape1 <= 0 || shape2 <= 0|| size <= 0) return R_NaN;
    /* FIXME: size = 0 is well defined ! */

    R_Q_P01_boundaries(p, 0, R_PosInf);

    double prob =  shape1 <= 1? R_PosInf : shape2 / (shape1 - 1); //TODO alternative method for shape1 <= 3
    Q = 1.0 / prob;
    P = (1.0 - prob) * Q;
    mu = size * prob;
    sigma = sqrt(size * (shape1 + size - 1) * shape2 * (shape1 + shape2 - 1) / ((shape1 - 2) * R_pow_di(shape1 - 1, 2))); //SD
    gamma = (shape1 + 2 * size - 1) * (shape1 + 2 * shape2 - 1) / ((shape1 - 3) * sqrt(size * (shape1 + size - 1) * shape2 * (shape1 + shape2 - 1) / (shape1 - 2))); //Skewsizeess

    /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower_tail || log_p) {
	p = R_DT_qIv(p); /* need check again (cancellation!): */
	if (p == R_DT_0) return 0;
	if (p == R_DT_1) return R_PosInf;
    }
    /* temporary hack --- FIXME --- */
    if (p + 1.01*DBL_EPSILON >= 1.) return R_PosInf;

    /* y := approx.value (Cornish-Fisher expansion) :  */
    z = qnorm(p, 0., 1., /*lower_tail*/TRUE, /*log_p*/FALSE);
    y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

    z = pbetanbinom(y, size, shape1, shape2, /*lower_tail*/TRUE, /*log_p*/FALSE);

    /* fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    /* If the C-F value is not too large a simple search is OK */
    if(y < 1e5) return do_search(y, &z, p, size, shape1, shape2, 1);
    /* Otherwise be a bit cleverer in the search */
    {
	double incr = floor(y * 0.001), oldincr;
	do {
	    oldincr = incr;
	    y = do_search(y, &z, p, size, shape1, shape2, incr);
	    incr = fmax2(1, floor(incr/100));
	} while(oldincr > 1 && incr > y*1e-15);
	return y;
    }
}

// double qbetanbinom_mu(double p, double size, double mu, int lower_tail, int log_p)
// {
// /* FIXME!  Implement properly!! (not losing accuracy for very large size (prob ~= 1)*/
//     return qbetanbinom(p, size, /* prob = */ size/(size+mu), lower_tail, log_p);
// }

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2006  The R Core Team
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
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rbetanbinom(double n, double p)
 *
 *  DESCRIPTION
 *
 *    Random variates from the negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York:Springer-Verlag. Page 480.
 *
 *  METHOD
 *
 *    Generate lambda as gamma with shape parameter n and scale
 *    parameter p/(1-p).  Return a Poisson deviate with mean lambda.
 */


double rbetanbinom(double size, double shape1, double shape2)
{
//     if(!R_FINITE(size) || size <= 0 || shape1 <= 0 || shape2 <= 0)
// 	return R_NaN;
    return rnbinom(size, rbeta(shape1, shape2));
}

// double rbetanbinom_mu(double size, double mu)
// {
//     if(!R_FINITE(size) || !R_FINITE(mu) || size <= 0 || mu < 0)
// 	ML_ERR_return_NAN;
//     return (mu == 0) ? 0 : rpois(rgamma(size, mu / size));
// }

/*
MeanBetaBinomialNegativaDist<- function(n, a, b){ # a i b són paràmetres Beta
  res<- n*b/(a-1)
  res[a <= 1]<- Inf
  
  return (res)
}

VarBetaBinomialNegativaDist<- function(n, a, b){ # a i b són paràmetres Beta
  res<- n*(a+n-1)*b*(a+b-1)/((a-2)*((a-1)^2))
  res[a <= 2]<- Inf
  
  return (res)
}

TrobaParamBetaNegBinomDist<- function(n, mu, sigma, DEBUG=FALSE){ # return(data.frame(a,b))
# Hi ha restriccions en l'espai mu ~ sigma: alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mu=n*b/(a-1) , sigma= n*(a+n-1)*b*(a+b-1)/((a-2)*((a-1)^2))], [a,b]);
  a<- (2 * n * sigma + mu * n^2 + (mu^2 - mu) * n - mu^2) / (n * sigma - mu * n - mu^2)
  b<- (mu * sigma + mu^2 * n + mu^3) / (n * sigma - mu * n - mu^2)
  a[which(a <= 2 | b <= 0)]<- NA
  b[which(is.na(a))]<- NA
  
  return (data.frame(a,b))
}
*/