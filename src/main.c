#include "probability.h"
#include <stdio.h>
#include <Rmath.h>
#include <R.h>
#include "dpq.h"

int main()
{
    double shape1, shape2, mean, var, prob, probLog;
    int n, lower_tail, log_p;
    
    shape1=1;
    shape2=2;
    n=100;
    prob=0.5;
    probLog= log(prob);
    
    mean= shape1 / (shape1 + shape2);
    var = shape1 * shape2 / (((shape1 + shape2)*(shape1 + shape2)) *(shape1+shape2+1.0));
//     char string[] = "mu= %f\tsigma= %f";
    lower_tail=TRUE;
    log_p=TRUE;

    printf("p= %f\tlog(p)= %f\tR_DT_qIv= %f", prob, probLog, R_DT_qIv(probLog));
    printf("mu= %f\tsigma= %f\tn=%i\n", mean, var, n);

    printf("dbb(2,%i,%.2f,%.2f)= %f\n", n, shape1, shape2, dbetabinom(2,n,shape1,shape2, FALSE));
    printf("dbb(2,%i,%.2f,%.2f)log= %f\n", n, shape1, shape2, exp(dbetabinom(2,n,shape1,shape2, TRUE)));
    
    printf("pbb(2,%i,%.2f,%.2f)= %f\n", n, shape1, shape2, pbetabinom(2,n,shape1,shape2, TRUE, FALSE));
    printf("pbb(2,%i,%.2f,%.2f)log= %f\n", n, shape1, shape2, exp(pbetabinom(2,n,shape1,shape2, TRUE, TRUE)));
    
    printf("qbb(.9,%i,%.2f,%.2f)= %.0f\n", n, shape1, shape2, qbetabinom(.9,n,shape1,shape2, TRUE, FALSE));
    printf("qbb(.9,%i,%.2f,%.2f)log= %.0f\n", n, shape1, shape2, qbetabinom(log(.9),n,shape1,shape2, TRUE, TRUE));
    
    return 0;
}
