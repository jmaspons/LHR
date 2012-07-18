/*
 * AUTHOR
 *  Copyright (C) 2012, Joan Maspons, joanmaspons@gmail.com
 *  based on negative binomial distribution code from The R Core Team
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
 *   Computes the negative binomial distribution. For integer size,
 *   this is probability of x failures before the nth success in a
 *   sequence of Bernoulli trials. We do not enforce integer size, since //TODO check statement
 *   the distribution is well defined for non-integers,
 *   and this can be useful for e.g. overdispersed discrete survival times.
 */
//TODO check alternative parameterization nbinom_mu

#include <R.h>
#include <Rmath.h>
#include "betaBinom.h"
#include "dpq.h"
#include "R_ext/Visibility.h"

double attribute_hidden dbetanbinom_raw(double x, double size, double shape1, double shape2, int give_log)
{
    double ans;
 
    if (x < 0 || !R_FINITE(x)) return R_D__0;
    
//     ans = pochhammer(shape1, size) * pochhammer(size, x) * pochhammer(shape2, x) / (gamma(x+1) * pochhammer(shape1 + shape2, size) * pochhammer(size + shape1 + shape2, x)); //factorial(x) = gamma(x+1)
    ans = gammafn(shape1 + size) / gammafn(shape1) * gammafn(size + x) / gammafn(size) * gammafn(shape2 + x) / gammafn(shape2) / (gammafn(x+1) * gammafn(shape1 + shape2 + size) / gammafn(shape1 + shape2) * gammafn(size + shape1 + shape2 + x) / gammafn(size + shape1 + shape2)); //factorial(x) = gamma(x+1)
if (!R_FINITE(ans)){
  REprintf("Numerator=%7g Denominator= %7g\n", gammafn(shape1 + size) / gammafn(shape1) * gammafn(size + x) / gammafn(size) * gammafn(shape2 + x) / gammafn(shape2), gammafn(x+1) * gammafn(shape1 + shape2 + size) / gammafn(shape1 + shape2) * gammafn(size + shape1 + shape2 + x) / gammafn(size + shape1 + shape2));
}
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

/*
 *  DESCRIPTION
 *
 *	The distribution function of the Beta negative binomial distribution.
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
    if (!lower_tail) ans = 1 - ans; //TODO check lower_tail -> R_D_Cval
    return ans;
}

/*
 *  DESCRIPTION
 *
 *	The quantile function of the Beta negative binomial distribution.
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

double qbetanbinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p)
{
    double pi, i = 0;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
	return p + size + shape1 + shape2;
#endif
    if (shape1 <= 0 || shape2 <= 0|| size < 0) return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    
    
    while (pi < p) {
	pi += dbetanbinom_raw(i, size, shape1, shape2, log_p);
	i++;
    }
    
    return i; //TODO check lower_tail -> R_D_Cval
}

/*
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