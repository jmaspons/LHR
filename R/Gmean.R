# Geometric mean
## TODO: check paper about geometric mean:
# Habib, E. A. (2012). Geometric mean for negative and zero values. International Journal of Research and Reviews in Applied Sciences, 11, 419-432.
#' Geometric mean
#'
#' @param gr 
#'
#' @return
#' @export
#'
#' @examples gr<- r(pop) | lambda(pop) # grow rate
Gmean<- function(gr){
  if (all(gr != 0)){
    res<- prod(gr)^(1/length(gr)) # buffer overflow for large numbers
    if (is.finite(res))  return (res)
    res<- exp(mean(log(gr))) # NaN for negative numbers
    if (!is.na(res)) return (res)
  }else{
    return (G(mean(gr), var(as.numeric(gr))))
  }
}


#' Geometric mean approximation from mean and variance
#'
#' @return a numeric vector with the geometric mean.
#' @references Stearns, S. C. (2000). Daniel Bernoulli (1738): evolution and economics under risk. Journal of Biosciences, 25(3), 221–228.
#' 
#' @export
G<- function(...){
  UseMethod("G")
}

#' @rdname G
#' 
#' @param mu mean.
#' @param sigma2 variance.
#' @examples G(mu=5, sigma2=.1)
#' @export
G.numeric<- function(mu, sigma2){
  return (mu - sigma2 / (2 * mu))
}

#' @rdname G
#' 
#' @param pop a \code{\link{discretePopSim}} object.
#' @param growRate a character string (\code{r} or \code{lambda}) indicating which growth rate to use.
#' @examples G(discre)
#' @export
G.discretePopSim<- function(pop, growRate=c("r", "lambda")[1], ...){
  Gres<- switch(growRate,
         r={rPop<- as.numeric(r(pop), ...); # avoid var() returns a var/cov matrix
            G(mean(rPop), var(rPop))},
         lambda={lambdaPop<- as.numeric(lambda(pop), ...);
            GLambda<- G(mean(lambdaPop), var(lambdaPop))})
  return (Gres)
} 
