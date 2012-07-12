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

#include "betaBinom.h"
#include "dpq.h"
#include <R.h>
#include <Rmath.h>

double attribute_hidden dbetabinom_raw(double x, double n, double a, double b, int give_log)
{
    double ans;

    if (x == 0 && n == 0) return R_D__1;
    if (x < 0 || x > n) return( R_D__0 );

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
