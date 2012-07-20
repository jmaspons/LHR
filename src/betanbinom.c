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
 *   sequence of Bernoulli trials. We do not enforce integer size, since
 *   the distribution is well defined for non-integers,
 *   and this can be useful for e.g. overdispersed discrete survival times.
 */

#include <R.h>
#include <Rmath.h>
#include "betaBinom.h"
#include "dpq.h"
#include "R_ext/Visibility.h"

double attribute_hidden dbetanbinom_raw(double x, double size, double shape1, double shape2, int give_log)
{
    if (x < 0 || !R_FINITE(x)) return R_D__0;

    double ans;
    if (log_p) {
        ans = lgammafn(shape1 + size) - lgammafn(shape1) + lgammafn(size + x) - lgammafn(size) + lgammafn(shape2 + x) - lgammafn(shape2)
              - (lgammafn(x+1) + lgammafn(shape1 + shape2 + size) - lgammafn(shape1 + shape2) + lgammafn(size + shape1 + shape2 + x) - lgammafn(size + shape1 + shape2));
    } else {
        //FIXME gammafn doesn't suport x > 171
        ans = gammafn(shape1 + size) / gammafn(shape1) * gammafn(size + x) / gammafn(size) * gammafn(shape2 + x) / gammafn(shape2) /
              (gammafn(x+1) * gammafn(shape1 + shape2 + size) / gammafn(shape1 + shape2) * gammafn(size + shape1 + shape2 + x) / gammafn(size + shape1 + shape2)); //factorial(x) = gamma(x+1)
    }

    return ans;
}

double dbetanbinom(double x, double size, double shape1, double shape2, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
        return x + size + shape1 + shape2;
#endif

    if (shape1 <= 0 || shape2 <= 0 || size < 0) return R_NaN;
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
    double ans = 0;
    for (int i = 0; i <= x; ++i) {
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

    return R_D_Lval(ans);
}

/*
 *  DESCRIPTION
 *
 *	The quantile function of the Beta negative binomial distribution.
 *
 *  NOTES
 *
 *	x = the number of failures before the n-th success
 */

double qbetanbinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(size) || ISNAN(shape1) || ISNAN(shape2))
        return p + size + shape1 + shape2;
#endif
    if (shape1 <= 0 || shape2 <= 0|| size < 0) return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);

    p = R_DT_val(p); //TODO check /* need check again (cancellation!): */

    double pi = 0, i = 0;

    while (pi < p) {
        pi += dbetanbinom_raw(i, size, shape1, shape2, log_p);
        i++;
    }

    return i;
}

/*
 *  DESCRIPTION
 *
 *    Random variates from the Beta negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 */


double rbetanbinom(double size, double shape1, double shape2)
{
    return rnbinom(size, rbeta(shape1, shape2));
}

