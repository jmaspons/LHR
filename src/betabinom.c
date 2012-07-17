/*
 * AUTHOR
 *  2012 Joan Maspons, joanmaspons@gmail.com
 *
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
 *   To compute the beta binomial probability, call dbetabinom(x,n, a, b).
 *   This checks for argument validity, and calls dbetabinom_raw().
 *
 *   dbetabinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbetabinom().
 *     (1) dbetabinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *     (2) Also does not check for 0 < a and 0 < b or NaN's.
 *         Do this in the calling function.
 */
//TODO unify var names a= shape1, b = shape2, n=size
#include <R.h>
#include <Rmath.h>
//#include "betaBinom.h"
#include "dpq.h"
#include "R_ext/Visibility.h"

double attribute_hidden dbetabinom_raw(double x, double n, double a, double b, int give_log)
{
    double ans;

    if (x == 0 && n == 0) return R_D__1;
    if (x < 0 || x > n) return R_D__0;

    ans = choose(n,x) * beta(x+a, n-x+b) / beta(a, b);

    return R_D_qIv(ans);
}

double dbetabinom(double x, double n, double a, double b, int give_log)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n) || ISNAN(a) || ISNAN(b)) return x + n + a + b;
#endif

    if (a <= 0 || b <= 0 || R_D_negInonint(n)) return R_NaN;
    R_D_nonint_check(x);
    if (x < 0 || !R_FINITE(x)) return R_D__0;

    n = R_D_forceint(n);
    x = R_D_forceint(x);
    return dbetabinom_raw(x, n, a, b, give_log);
}

/*
 *
 *  DESCRIPTION
 *
 *    The distribution function of the beta binomial distribution.
 */

double attribute_hidden pbetabinom_raw(double x, double n, double a, double b, int log_p)
{
    int i;
    double ans = 0;
    for (i = 0; i <= x; ++i) {
        ans += dbetabinom_raw(i, n, a, b, log_p);
    }
    return ans;
}

double pbetabinom(double x, double n, double a, double b, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(a) || ISNAN(b))
        return x + n + a + b;
    if (!R_FINITE(n)) return R_NaN; //TODO || !R_FINITE(shape1) || !R_FINITE(shape2)
#endif

    if (a <= 0 || b <= 0 || R_D_nonint(n)) return R_NaN;
    n = R_D_forceint(n);

    if (n < 0) return R_NaN;

    if (x < 0) return R_DT_0;
    x = floor(x + 1e-7);
    if (n <= x) return R_DT_1;

    double ans = pbetabinom_raw(x, n, a, b, log_p);
    if (!lower_tail) ans = 1 - ans;
    return ans;
}
/*
 *
 *  DESCRIPTION
 *
 *	The quantile function of the Beta binomial distribution.
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation. A search is then conducted of values close to
 *	this initial start point.
 */

static double
do_search(double y, double *z, double p, double n, double a, double b, double incr)
{
    if(*z >= p) {
        /* search to the left */
#if DEBUG_c_qbetabinom == 1
        printf("\tnew z=%7g >= p = %7g  --> search to left (y--) ..\n", *z,p);
#endif
#ifdef DEBUG_qbetabinom
        REprintf("\tnew z=%7g >= p = %7g  --> search to left (y--) ..\n", z,p);
#endif
        for(;;) {
            double newz;
            if(y == 0 ||
                    (newz = pbetabinom(y - incr, n, a, b, /*l._t.*/TRUE, /*log_p*/FALSE)) < p)
                return y;
            y = fmax2(0, y - incr);
            *z = newz;
#if DEBUG_c_qbetabinom == 1
	    printf("\tnew z=%7g >= p = %7g  --> new y=%.0f ..\n", *z,p, y);
#endif    
        }
    }
    else {		/* search to the right */
#if DEBUG_c_qbetabinom == 1
        printf("\tnew z=%7g < p = %7g  --> search to right (y++) ..\n", *z,p);
#endif
#ifdef DEBUG_qbetabinom
        REprintf("\tnew z=%7g < p = %7g  --> search to right (y++) ..\n", z,p);
#endif
        for(;;) {
            y = fmin2(y + incr, n);
            if(y == n ||
                    (*z = pbetabinom(y, n, a, b, /*l._t.*/TRUE, /*log_p*/FALSE)) >= p)
                return y;
#if DEBUG_c_qbetabinom == 1
	    printf("\tnew z=%7g >= p = %7g  --> new y=%.0f ..\n", *z,p, y);
#endif 
        }
    }
}

double qbetabinom(double p, double n, double a, double b, int lower_tail, int log_p)
{
    double mu, pr, sigma, gamma, z, y;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(n) || ISNAN(a) || ISNAN(b))
        return p + n + a + b;
#endif
    if(!R_FINITE(n)) return R_NaN;
    /* if log_p is true, p = -Inf is a legitimate value */
    if(!R_FINITE(p) && !log_p) return R_NaN;

    if(n != floor(n + 0.5)) return R_NaN;
    if (a <= 0 || b <= 0 || n < 0) return R_NaN;

    R_Q_P01_boundaries(p, 0, n);

    pr = a / (a + b); // Beta Mean
    mu = n * pr; // Mean
    sigma = sqrt(n * a * b * (n + a + b) / (a + b + 1)) / (b + a); // SD
    gamma = (a + b + 2 * n) * (b - a) / (a + b + 2) * sqrt((1 + a + b) / (n * a * b * (n + a + b))); // Skewness

#if DEBUG_c_qbetabinom == 1
    printf("qbetabinom(p=%7g, n=%g, pr=%7g, l.t.=%d, log=%d): sigm=%g, gam=%g\n",
             p,n,pr, lower_tail, log_p, sigma, gamma);
#endif   
#ifdef DEBUG_qbetabinom
    REprintf("qbetabinom(p=%7g, n=%g, pr=%7g, l.t.=%d, log=%d): sigm=%g, gam=%g\n",
             p,n,pr, lower_tail, log_p, sigma, gamma);
#endif
    /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower_tail || log_p) {
        p = R_DT_qIv(p); /* need check again (cancellation!): */
        if (p == 0.) return 0.;
        if (p == 1.) return n;
    }
    /* temporary hack --- FIXME --- */
    if (p + 1.01*DBL_EPSILON >= 1.) return n;

    /* y := approx.value (Cornish-Fisher expansion) :  */
    z = qnorm(p, 0., 1., /*lower_tail*/TRUE, /*log_p*/FALSE);
    y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

    if(y > n) /* way off */ y = n;
#if DEBUG_c_qbetabinom == 1
    printf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1-p, z, y);
#endif
#ifdef DEBUG_qbetabinom
    REprintf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1-p, z, y);
#endif
    z = pbetabinom(y, n, a, b, /*lower_tail*/TRUE, /*log_p*/FALSE);

    /* fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    if(n < 1e5) return do_search(y, &z, p, n, a, b, 1);
    /* Otherwise be a bit cleverer in the search */
    {
        double incr = floor(n * 0.001), oldincr;
        do {
            oldincr = incr;
            y = do_search(y, &z, p, n, a, b, incr);
            incr = fmax2(1, floor(incr/100));
        } while(oldincr > 1 && incr > n*1e-15);
        return y;
    }
}


double rbetabinom(double n, double a, double b)
{
    return rbinom(n, rbeta(a, b));
}
