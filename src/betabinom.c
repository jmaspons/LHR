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

double attribute_hidden dbetabinom_raw(double x, double size, double shape1, double shape2, int give_log)
{
    double ans;

    if (x == 0 && size == 0) return R_D__1;
    if (x < 0 || x > size) return R_D__0;

    ans = choose(size,x) * beta(x+shape1, size-x+shape2) / beta(shape1, shape2);

    return R_D_qIv(ans);
}

double dbetabinom(double x, double size, double shape1, double shape2, int give_log)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2)) return x + size + shape1 + shape2;
#endif

    if (shape1 <= 0 || shape2 <= 0 || R_D_negInonint(size)) return R_NaN;
    R_D_nonint_check(x);
    if (x < 0 || !R_FINITE(x)) return R_D__0;

    size = R_D_forceint(size);
    x = R_D_forceint(x);
    return dbetabinom_raw(x, size, shape1, shape2, give_log);
}

/*
 *
 *  DESCRIPTION
 *
 *    The distribution function of the beta binomial distribution.
 */

double attribute_hidden pbetabinom_raw(double x, double size, double shape1, double shape2, int log_p)
{
    int i;
    double ans = 0;
    for (i = 0; i <= x; ++i) {
        ans += dbetabinom_raw(i, size, shape1, shape2, log_p);
    }
    return ans;
}

double pbetabinom(double x, double size, double shape1, double shape2, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
        return x + size + shape1 + shape2;
    if (!R_FINITE(size)) return R_NaN; //TODO || !R_FINITE(shape1) || !R_FINITE(shape2)
#endif

    if (shape1 <= 0 || shape2 <= 0 || R_D_nonint(size)) return R_NaN;
    size = R_D_forceint(size);

    if (size < 0) return R_NaN;

    if (x < 0) return R_DT_0;
    x = floor(x + 1e-7);
    if (size <= x) return R_DT_1;

    double ans = pbetabinom_raw(x, size, shape1, shape2, log_p);
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
do_search(double y, double *z, double p, double size, double shape1, double shape2, double incr)
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
                    (newz = pbetabinom(y - incr, size, shape1, shape2, /*l._t.*/TRUE, /*log_p*/FALSE)) < p)
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
            y = fmin2(y + incr, size);
            if(y == size ||
                    (*z = pbetabinom(y, size, shape1, shape2, /*l._t.*/TRUE, /*log_p*/FALSE)) >= p)
                return y;
#if DEBUG_c_qbetabinom == 1
	    printf("\tnew z=%7g >= p = %7g  --> new y=%.0f ..\n", *z,p, y);
#endif 
        }
    }
}

double qbetabinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p)
{
    double mu, pr, sigma, gamma, z, y;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
        return p + size + shape1 + shape2;
#endif
    if(!R_FINITE(size)) return R_NaN;
    /* if log_p is true, p = -Inf is a legitimate value */
    if(!R_FINITE(p) && !log_p) return R_NaN;

    if(size != floor(size + 0.5)) return R_NaN;
    if (shape1 <= 0 || shape2 <= 0 || size < 0) return R_NaN;

    R_Q_P01_boundaries(p, 0, size);

    pr = shape1 / (shape1 + shape2); // Beta Mean
    mu = size * pr; // Mean
    sigma = sqrt(size * shape1 * shape2 * (size + shape1 + shape2) / (shape1 + shape2 + 1)) / (shape2 + shape1); // SD
    gamma = (shape1 + shape2 + 2 * size) * (shape2 - shape1) / (shape1 + shape2 + 2) * sqrt((1 + shape1 + shape2) / (size * shape1 * shape2 * (size + shape1 + shape2))); // Skewness

#if DEBUG_c_qbetabinom == 1
    printf("qbetabinom(p=%7g, size=%g, pr=%7g, l.t.=%d, log=%d): sigm=%g, gam=%g\n",
             p,size,pr, lower_tail, log_p, sigma, gamma);
#endif   
#ifdef DEBUG_qbetabinom
    REprintf("qbetabinom(p=%7g, size=%g, pr=%7g, l.t.=%d, log=%d): sigm=%g, gam=%g\n",
             p,size,pr, lower_tail, log_p, sigma, gamma);
#endif
    /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower_tail || log_p) {
        p = R_DT_qIv(p); /* need check again (cancellation!): */
        if (p == 0.) return 0.;
        if (p == 1.) return size;
    }
    /* temporary hack --- FIXME --- */
    if (p + 1.01*DBL_EPSILON >= 1.) return size;

    /* y := approx.value (Cornish-Fisher expansion) :  */
    z = qnorm(p, 0., 1., /*lower_tail*/TRUE, /*log_p*/FALSE);
    y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

    if(y > size) /* way off */ y = size;
#if DEBUG_c_qbetabinom == 1
    printf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1-p, z, y);
#endif
#ifdef DEBUG_qbetabinom
    REprintf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1-p, z, y);
#endif
    z = pbetabinom(y, size, shape1, shape2, /*lower_tail*/TRUE, /*log_p*/FALSE);

    /* fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    if(size < 1e5) return do_search(y, &z, p, size, shape1, shape2, 1);
    /* Otherwise be a bit cleverer in the search */
    {
        double incr = floor(size * 0.001), oldincr;
        do {
            oldincr = incr;
            y = do_search(y, &z, p, size, shape1, shape2, incr);
            incr = fmax2(1, floor(incr/100));
        } while(oldincr > 1 && incr > size*1e-15);
        return y;
    }
}


double rbetabinom(double size, double shape1, double shape2)
{
    return rbinom(size, rbeta(shape1, shape2));
}
