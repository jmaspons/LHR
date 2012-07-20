/*
 * AUTHOR
 *  Copyright (C) 2012, Joan Maspons, joanmaspons@gmail.com
 *  based on binomial distribution code from The R Core Team
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
 *   To compute the beta binomial probability, call dbetabinom(x, size, shape1, shape2).
 *   This checks for argument validity, and calls dbetabinom_raw().
 *
 *   dbetabinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbetabinom().
 *     (1) dbetabinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *     (2) Also does not check for 0 < shape1 and 0 < shape2 or NaN's.
 *         Do this in the calling function.
 */

#include <R.h>
#include <Rmath.h>
//#include "betaBinom.h"
#include "dpq.h"
#include "R_ext/Visibility.h"

double attribute_hidden dbetabinom_raw(double x, double size, double shape1, double shape2, int give_log)
{
    if (x == 0 && size == 0) return R_D__1;
    if (x < 0 || x > size) return R_D__0;

    if (log_p) {
        return (lchoose(size,x) + lbeta(x+shape1, size-x+shape2) - lbeta(shape1, shape2));
    } else {
        //FIXME for choose(x,y) => max x = 1030 if y = x/2 & x>y
        //FIXME for beta(x,y) => max(x,y) = 509
        return (choose(size,x) * beta(x+shape1, size-x+shape2) / beta(shape1, shape2));
    }
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
 *    The distribution function of the Beta binomial distribution.
 */

double attribute_hidden pbetabinom_raw(double x, double size, double shape1, double shape2, int log_p)
{
    double ans = 0;
    for (int i = 0; i <= x; ++i) {
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

    return R_D_Lval(ans);
}
/*
 *
 *  DESCRIPTION
 *
 *	The quantile function of the Beta binomial distribution.
 */

double qbetabinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p)
{
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

    p = R_DT_val(p);//TODO check /* need check again (cancellation!): */

    double pi = 0;

    for (int i = 0; i < size; ++i) {
        pi += dbetabinom_raw(i, size, shape1, shape2, log_p);
        if (pi > p) return i;
    }
    return size;
}


double rbetabinom(double size, double shape1, double shape2)
{
    return rbinom(size, rbeta(shape1, shape2));
}
