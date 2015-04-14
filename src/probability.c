#include <Rmath.h>

int sdistri(double * Pdist, double * Vdist, int n, double replace, double ** res){
    int i, found;
    double * resPtr, cumsum;
    
    for (i = 0, found = 0, cumsum = 0; i < n; i++) {
      resPtr[0] += Pdist[i] * Vdist[i]; // mean = sum(P(i) * i / n0))
      if (!found) {
        cumsum += Pdist[i];
        if (Vdist[i] >= replace) {
          resPtr[3] = 1 - cumsum; // pReplace = P(i) | i == 2 * n0
          found = 1;
        }
      }
    }

    for (i = 0; i <= n; i++) {
      resPtr[1] += Pdist[i] * R_pow_di(Vdist[i] - resPtr[0], 2); // var = sum(P(i) * (i / n0 - mean)^2)
    }

    resPtr[2] = resPtr[0] - 2 * resPtr[1] / resPtr[0]; // G = mean - 2 * var / mean
    
    * res = resPtr;

    return 4;
}
