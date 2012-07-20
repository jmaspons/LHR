#include "betaBinom.h"
#include <stdio.h>
#include <Rmath.h>
#include <R.h>

int main()
{
    double shape1, shape2, mean, var;
    int n;
    shape1=1;
    shape2=2;
    n=100;

    mean= shape1 / (shape1 + shape2);
    var = shape1 * shape2 / (((shape1 + shape2)*(shape1 + shape2)) *(shape1+shape2+1.0));
    char string[] = "mu= %f\tsigma= %f";

    printf("mu= %f\tsigma= %f\tn=%i\n", mean, var, n);

    printf("dbb(2,%i,%.2f,%.2f)= %f\n", n, shape1, shape2, dbetabinom(2,n,shape1,shape2, FALSE));
//   GetRNGstate();
    printf("rbb(%i,%.2f,%.2f)= %.0f\n", n, shape1, shape2, rbetabinom(n,shape1,shape2));
//   PutRNGstate();
    printf("pbb(2,%i,%.2f,%.2f)= %f\n", n, shape1, shape2, pbetabinom(2,n,shape1,shape2, TRUE, FALSE));
    printf("qbb(.9,%i,%.2f,%.2f)= %.0f\n", n, shape1, shape2, qbetabinom(.9,n,shape1,shape2, TRUE, FALSE));

    return 0;
}
