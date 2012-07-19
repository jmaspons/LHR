
dbetabinom <- function(x, size, shape1, shape2, mean, variance, log = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "dbetabinom", x, size, betaPar[[1]], betaPar[[2]], log)
    }
    else .External("actuar_do_dpq", "dbetabinom", x, size, shape1, shape2, log)
}

pbetabinom <- function(q, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "pbetabinom", q, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "pbetabinom", q, size, shape1, shape2, lower.tail, log.p)
}

qbetabinom <- function(p, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "qbetabinom", p, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "qbetabinom", p, size, shape1, shape2, lower.tail, log.p)
}

rbetabinom <- function(n, size, shape1, shape2, mean, variance)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaBinomDist(size, mean, variance)
        .External("actuar_do_random", "rbetabinom", n, size, betaPar[[1]], betaPar[[2]])
    }
    else .External("actuar_do_random", "rbetabinom", n, size, shape1, shape2)
}

MeanBetaBinomialDist<- function(n, a, b){ # a i b són paràmetres Beta
  return (n*a/(a+b))
}

VarBetaBinomialDist<- function(n, a, b){ # a i b són paràmetres Beta
  return (n*a*b*(a+b+n)/(((a+b)^2)*(a+b+1)))
}

TrobaParamBetaBinomDist<- function(n, mu, sigma){ # return(data.frame(a,b))
# Hi ha restriccions en l'espai mu ~ sigma: alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mu=n*a/(a+b) , sigma= n*a*b*(a+b+n)/(((a+b)^2)*(a+b+1))], [a,b]);
  a<- -(mu * sigma - mu^2 * n + mu^3) / (n * sigma - mu * n + mu^2)
  b<- -((n - mu) * sigma - mu * n^2 + 2 * mu^2 * n - mu^3)/ (n * sigma - mu * n + mu^2)
  a[which(a <= 0 | b <= 0)]<- NA
  b[which(is.na(a))]<- NA
  
  return (data.frame(a,b))
}