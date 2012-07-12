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
 *  DESCRIPTION
 *
 *    The distribution function of the beta binomial distribution.
 */
#include "betaBinom.h"
#include "dpq.h"
#include <R.h>
#include <Rmath.h>

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
    if (!R_FINITE(n)) return R_NaN;
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
