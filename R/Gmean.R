# Geometric mean
## TODO: check paper about geometric mean:
# Habib, E. A. (2012). Geometric mean for negative and zero values. International Journal of Research and Reviews in Applied Sciences, 11, 419-432.
#' Geometric mean
#'
#' @param x 
#' @param na.rm
#'
#' @return the geometric mean
#' @export
Gmean<- function(x, na.rm=TRUE, ...) UseMethod("Gmean")

#' @rdname Gmean
#' @examples
#' pop<- discretePopSim(b=1, j=.5, a=.5)
#' Gmean(pop)
#' @export
Gmean.discretePopSim<- function(x, na.rm=TRUE, ...){
  # NextMethod(object=lambda(x))
  Gmean(lambda(x))
}

Gmean.matrix<- function(x, na.rm=TRUE, ...){
  if (all(x != 0)){
    res<- prod(x)^(1/length(x)) # buffer overflow for large numbers
    if (is.finite(res))  return (res)
    res<- exp(mean(log(x))) # NaN for negative numbers
    if (!is.na(res)) return (res)
  }else{
    return (G(mean(x, na.rm=na.rm), var(as.numeric(x), na.rm=na.rm)))
  }
}


#' Geometric mean approximation from mean and variance
#'
#' @return a numeric vector with the geometric mean.
#' @references Stearns, S. C. (2000). Daniel Bernoulli (1738): evolution and economics under risk. Journal of Biosciences, 25(3), 221-228.
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
G.numeric<- function(mu, sigma2, ...){
  return (mu - sigma2 / (2 * mu))
}

#' @rdname G
#' 
#' @param pop a \code{\link{discretePopSim}} object.
#' @param growRate a character string (\code{r} or \code{lambda}) indicating which growth rate to use.
#' @examples
#' pop<- discretePopSim(b=1, j=.5, a=.5)
#' G(pop)
#' @export
G.discretePopSim<- function(pop, growRate=c("lambda", "r"), na.rm=TRUE, ...){
  growRate<- match.arg(growRate)
  Gres<- switch(growRate,
         r={rPop<- as.numeric(as.matrix(r(pop))); # avoid var() returns a var/cov matrix
            G(mean(rPop, na.rm=na.rm), var(rPop, na.rm=na.rm))},
         lambda={lambdaPop<- as.numeric(as.matrix(lambda(pop)));
            GLambda<- G(mean(lambdaPop, na.rm=na.rm), var(lambdaPop, na.rm=na.rm))})
  return (Gres)
} 
