#include <R.h>
#include <Rmath.h>

void fbeta(double mean, double var, double shape[2])
{

#ifdef IEEE_754
    if (ISNAN(mean) || ISNAN(var)) {
        shape[0] = R_NaN;
        shape[1] = R_NaN;
        return;
    }
#endif

    shape[0] = -(mean * var + R_pow_di(mean, 3) - R_pow_di(mean, 2)) / var;
    shape[1] = ((mean-1) * var + R_pow_di(mean, 3) - 2 * R_pow_di(mean, 2) + mean) / var;

    if (shape[0] <= 0 || shape[1] <= 0) {
        shape[0] = R_NaN;
        shape[1] = R_NaN;
    }

    return;
}
