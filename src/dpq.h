/* Utilities for `dpq' handling (density/probability/quantile) from actuar package*/

/* give_log in "d" & "mgf";  log_p in "p" & "q" : */
#define give_log log_p

#define R_D__0  (log_p ? R_NegInf : 0.)
#define R_D__1  (log_p ? 0. : 1.)
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)

/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
#define R_D_Lval(p)     (lower_tail ? (p) : (0.5 - (p) + 0.5))  /*  p  */
#define R_D_Cval(p)     (lower_tail ? (0.5 - (p) + 0.5) : (p))  /*  1 - p */

#define R_D_val(x)      (log_p  ? log(x) : (x))         /*  x  in pF(x,..) */
#define R_D_qIv(p)      (log_p  ? exp(p) : (p))         /*  p  in qF(p,..) */
#define R_D_exp(x)      (log_p  ?  (x)   : exp(x))      /* exp(x) */
#define R_D_log(p)      (log_p  ?  (p)   : log(p))      /* log(p) */
#define R_D_Clog(p)     (log_p  ? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

#define R_DT_val(x)     (lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)    (lower_tail ? R_D_Clog(x) : R_D_val(x))
//#define R_DT_qIv(p)      R_D_Lval(R_D_qIv(p))		 /*  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))
/*Boundaries*/
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)          \
    if (log_p) {                                        \
        if(p > 0)                                       \
            return R_NaN;				\
        if(p == 0) /* upper bound*/                     \
            return lower_tail ? _RIGHT_ : _LEFT_;       \
        if(p == R_NegInf)                               \
            return lower_tail ? _LEFT_ : _RIGHT_;       \
    }                                                   \
    else { /* !log_p */                                 \
        if(p < 0 || p > 1)                              \
            return R_NaN;				\
        if(p == 0)                                      \
            return lower_tail ? _LEFT_ : _RIGHT_;       \
        if(p == 1)                                      \
            return lower_tail ? _RIGHT_ : _LEFT_;       \
    }

/* Infinite limit in "lev" */
#define R_VG__0(x, y)   (R_FINITE(x) ? R_pow(x, y) : 0.)

#define R_NaN  (0.0/0.0)


/* Code from R core dpq.h */

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

/* Rprintf instead of MATHLIB_WARNING from "nmath.h"*/
#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	Rprintf("non-integer x = %f", x);	\
	return R_D__0;					\
   }

/* Extras */
#define pochhammer(x, n)	gammafn(x+n)/gammafn(x)
