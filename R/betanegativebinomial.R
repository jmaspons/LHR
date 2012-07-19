
dbetanbinom <- function(x, size, shape1, shape2, mean, variance, log = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaNegBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "dbetanbinom", x, size, betaPar[[1]], betaPar[[2]], log)
    }
    else .External("actuar_do_dpq", "dbetanbinom", x, size, shape1, shape2, log)
}

pbetanbinom <- function(q, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaNegBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "pbetanbinom", q, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "pbetanbinom", q, size, shape1, shape2, lower.tail, log.p)
}

qbetanbinom <- function(p, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaNegBinomDist(size, mean, variance)
        .External("actuar_do_dpq", "qbetanbinom", p, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "qbetanbinom", p, size, shape1, shape2, lower.tail, log.p)
}

rbetanbinom <- function(n, size, shape1, shape2, mean, variance)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- TrobaParamBetaNegBinomDist(size, mean, variance)
        .External("actuar_do_random", "rbetanbinom", n, size, betaPar[[1]], betaPar[[2]])
    }
    else .External("actuar_do_random", "rbetanbinom", n, size, shape1, shape2)
}

MeanBetaBinomialNegativaDist<- function(n, a, b){ # a i b són paràmetres Beta
  res<- n*b/(a-1)
  res[a <= 1]<- Inf
  
  return (res)
}

VarBetaBinomialNegativaDist<- function(n, a, b){ # a i b són paràmetres Beta
  res<- n*(a+n-1)*b*(a+b-1)/((a-2)*((a-1)^2))
  res[a <= 2]<- Inf
  
  return (res)
}

TrobaParamBetaNegBinomDist<- function(n, mu, sigma){ # return(data.frame(a,b))
# Hi ha restriccions en l'espai mu ~ sigma: alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mu=n*b/(a-1) , sigma= n*(a+n-1)*b*(a+b-1)/((a-2)*((a-1)^2))], [a,b]);
  a<- (2 * n * sigma + mu * n^2 + (mu^2 - mu) * n - mu^2) / (n * sigma - mu * n - mu^2)
  b<- (mu * sigma + mu^2 * n + mu^3) / (n * sigma - mu * n - mu^2)
  a[which(a <= 2 | b <= 0)]<- NA
  b[which(is.na(a))]<- NA
  
  return (data.frame(a,b))
}