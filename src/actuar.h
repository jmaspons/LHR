#include <Rinternals.h>

/*Error messages */
#define R_MSG_NA        _("NaNs produced")

/* Functions accessed from .External() */
SEXP actuar_do_dpq(SEXP args);
SEXP actuar_do_dpq1(int code, SEXP args);
SEXP actuar_do_dpq2(int code, SEXP args);
SEXP actuar_do_dpq3(int code, SEXP args);
SEXP actuar_do_dpq4(int code, SEXP args);
SEXP actuar_do_dpq5(int code, SEXP args);

SEXP actuar_do_random(SEXP args);
SEXP actuar_do_random1(int code, SEXP args);
SEXP actuar_do_random2(int code, SEXP args);
SEXP actuar_do_random3(int code, SEXP args);
SEXP actuar_do_random4(int code, SEXP args);

/* Utility functions */
/*   Matrix algebra */
void actuar_expm(double *x, int n, double *z);
double actuar_expmprod(double *x, double *M, double *y, int n);
void actuar_matpow(double *x, int n, int k, double *z);
void actuar_solve(double *A, double *B, int n, int p, double *z);

/*   Sampling */
int SampleSingleValue(int n, double *p);

/*   Three parameter distributions, hence associated with dpq3 */
double dbetabinom(double x, double size, double shape1, double shape2, int give_log);
double pbetabinom(double q, double size, double shape1, double shape2, int lower_tail, int log_p);
double qbetabinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p);
double rbetabinom(double size, double shape1, double shape2);
// double mbetabinom(double order, double size, double shape1, double shape2, int give_log);
// double levbetabinom(double limit, double size, double shape1, double shape2, double order, int give_log);

// double dbetanbinom(double x, double size, double shape1, double shape2, int give_log);
// double pbetanbinom(double q, double size, double shape1, double shape2, int lower_tail, int log_p);
// double qbetanbinom(double p, double size, double shape1, double shape2, int lower_tail, int log_p);
// double rbetanbinom(double size, double shape1, double shape2);
// double mbetanbinom(double order, double size, double shape1, double shape2, int give_log);
// double levbetanbinom(double limit, double size, double shape1, double shape2, double order, int give_log);

/* Definitions for the table linking the first group of functions to
 * the second one. Table found in names.c */
typedef struct {
    char *name;
    SEXP (*cfun)(int, SEXP);
    int code;
} FUNTAB;
extern FUNTAB fun_tab[];
