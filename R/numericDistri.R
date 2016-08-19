# Functions for compound discrete probability distributions resolved numerically
# numericDistri class

#' Numeric representation of discrete distributions
#' 
#' @name numericDistri
#' @exportClass numericDistri
NULL

## TODO: call distriBeta* or distri* according to p parameter (p = c(shape1, shape2) | p)

## Binomial distribution ----
#' @rdname numericDistri
#'
#' @inheritParams distriBinom
#' @param size 
#' @param prob 
#' @param logP 
#'
#' @return
#' @importFrom stats dbinom dnbinom qbinom qnbinom
#' @export
distriBinom<- function(size, prob, logP=FALSE){
  UseMethod("distriBinom")
}

#' @rdname numericDistri
#' @export
distriBinom.numeric<- function(size, prob, logP=FALSE){
  res<- data.frame(x=0:size, p=dbinom(x=0:size, size, prob, log=logP))
  attributes(res)$p.omitted<- 0
  attributes(res)$parameters<- list(size=size, prob=prob)
  attributes(res)$logP<- logP
  class(res)<- c("binom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# size: numericDistri object which compound a binomial distribution as a size parameter.
# p: probability
#' @rdname numericDistri
#' @export
distriBinom.numericDistri<- function(size, prob, logP=FALSE){
  if (logP != attributes(size)$logP){ # Transform p to log scale if necessary
    size<- logP(size, logP=logP)
  }
  
  res<- .External("binomialCompound", size$x, size$p, prob, logP, max(size$x))
  res<- data.frame(x=0:(length(res) - 1), p=res)
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(size=c(list(distribution=class(size)[1]), attributes(size)$parameters), prob=prob)
  attributes(res)$logP<- logP
  class(res)<- "compoundBinom"
  if (inherits(size, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}


## Beta binomial distribution ----
#' @rdname numericDistri
#'
#' @inheritParams distriBinom
#' @param shape1 
#' @param shape2 
#'
#' @return
#' @export
#'
#' @examples
distriBetaBinom<- function(size, shape1, shape2, logP=FALSE){
  UseMethod("distriBetaBinom")
}

#' @rdname numericDistri
#' @export
distriBetaBinom.numeric<- function(size, shape1, shape2, logP=FALSE){
  res<- data.frame(x=0:size, p=dbetabinom(x=0:size, size, shape1, shape2, log=logP))
  attributes(res)$p.omitted<- 0
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var)
  attributes(res)$logP<- logP
  class(res)<- c("betaBinom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# size: numericDistribution object which compound a binomial distribution as a size parameter.
# p: probability
#' @rdname numericDistri
#' @export
distriBetaBinom.numericDistri<- function(size, shape1, shape2, logP=FALSE){
  if (logP != attributes(size)$logP){ # Transform p to log scale if necessary
    size<- logP(size, logP=logP)
  }
  
  res<- .External("BetaBinomialCompound", size$x, size$p, shape1, shape2, logP, max(size$x))
  res<- data.frame(x=0:(length(res) - 1), p=res)
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=c(list(distribution=class(size)[1]), attributes(size)$parameters), prob=prob$mean, var=prob$var)
  attributes(res)$logP<- logP
  class(res)<- "compoundBetaBinom"
  if (inherits(size, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}

## Negative binomial distribution ----
#' @rdname numericDistri
#' 
#' @inheritParams distriBinom
#' @param mu 
#' @param maxPomitted
#' @param maxX 
#' @return
#' @export
#'
#' @examples
distriNegBinom<- function(size, prob, mu, logP=FALSE, maxPomitted=0.01, maxX){
  UseMethod("distriNegBinom")
}

#' @rdname numericDistri
#' @export  
distriNegBinom.numeric<- function(size, prob, mu, logP=FALSE, maxPomitted=0.01, maxX){
  if (missing(prob) & !missing(mu)) prob<- size / (size+mu) #alternative parametrization ?dnbinom
  if (missing(maxX)){
    maxX<- qnbinom(maxPomitted, size=size, prob=prob, lower.tail=FALSE, log.p=logP)
  }
  res<- data.frame(x=0:maxX, p=dnbinom(x=0:maxX, size=size, prob=prob, log=logP))
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(size=size, prob=prob)
  attributes(res)$logP<- logP
  class(res)<- c("nbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}

#TODO distriNegBinom.numericDistri

## Beta negative binomial distribution ----
#' @rdname numericDistri
#' @export
distriBetaNegBinom<- function(size, shape1, shape2, logP=FALSE, maxPomitted=0.01, maxX){
  if (missing(maxX)){
    maxX<- qbetanbinom(maxPomitted, size=size, shape1=shape1, shape2=shape2, lower.tail=FALSE, log.p=logP)
  }
  res<- data.frame(x=0:maxX, p=dbetanbinom(x=0:maxX, size=size, shape1=shape1, shape2=shape2, log=logP))
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var)
  attributes(res)$logP<- logP
  class(res)<- c("betaNbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}

#TODO distriBetaNegBinom.numericDistri

## Change probability scale (p vs logp) ----
#' @rdname numericDistri
#' @export
logP<- function(distri, logP=TRUE){
  if (attributes(distri)$logP == logP) return (distri)
  
  if (logP){
    distri$p<- log(distri$p)
  }else{
    distri$p<- exp(distri$p)
  }
  attributes(distri)$logP<- logP
  return (distri)
}


## Generic methods for numericDistri class ----

#' @export
print<- function(x, ...) UseMethod("print")

#' @export
print.numericDistri<- function(x, ...){
  cat("\t", class(x)[1], "distribution\n")
  cat("Probability omitted: ", attributes(x)$p.omitted, "\n")
  if (attributes(x)$logP) cat("p in Log probability scale\n")
  cat("Parameters:\n")
  utils::str(attributes(x)$parameters)
  cat("\n")
  print(utils::head(data.frame(x), n=15), ...)
  if (nrow(x) > 15) cat("\t ...\t", nrow(x) - 15, "rows omited.\n")
}

#' @rdname numericDistri
#' @export
summary.numericDistri<- function(object, ...){
  cat("\t", class(object)[1], "distribution\n")
  cat("Probability omitted: ", attributes(object)$p.omitted, "\nParameters:\n")
  if (attributes(object)$logP) cat("p in Log probability scale\n")
  utils::str(attributes(object)$parameters)
  cat("\n")
  res<- sdistri(object)
  print(res)
  invisible(res)
}

#' @rdname numericDistri
#' @export
plot.numericDistri<- function(x, y, cum=FALSE, ...){
  if (cum){
    distri<- cumsum(x)
    distri$p<- distri$cump
  }
  type<- ifelse(cum, "s", "p")
  graphics::plot.default(distri, type=type, ...)
}


## Stats ----
#' @rdname numericDistri
#' @export
mean.numericDistri<- function(x, ...){
  if(attributes(x)$logP){
      x$p<- exp(x$p)
  }
  
  res<- stats::weighted.mean(x$x, x$p, ...)

  return (res)
}


#' @export
var<- function(x, ...){
  UseMethod("var")
}

#' @export
var.default<- function(x, ...) stats::var(x, ...)

#' @rdname numericDistri
#' @export
var.numericDistri<- function(x, ...){
  distri<- logP(x, logP=FALSE)
  sum(distri$p * (distri$x - mean(distri))^2)
}

#' @export
sdistri<- function(distri){
  UseMethod("sdistri")
}

#' @rdname numericDistri
#' @export
sdistri.numericDistri<- function(distri){
  distri<- logP(distri, logP=FALSE)

  meanD<- mean(distri)
  varD<- sum(distri$p * (distri$x - meanD)^2) #var(distri)
  GmeanD<- G(meanD, varD)
  
  return (data.frame(Gmean=GmeanD, mean=meanD, var=varD))
}


#' @export
cumsum<- function(x, ...){
  UseMethod("cumsum")
}

#' @export
cumsum.default<- function(x, ...) base::cumsum(x)

#' @rdname numericDistri
#' @export
cumsum.numericDistri<- function(x, ...){
  if(attributes(x)$logP){
    p<- exp(x$p)
    x$cump<- cumsum(p)
  }else{
    x$cump<- cumsum(x$p)
  }
  
  return(x)
}


#' Random deviates from a numericDistribution
#' 
#' @export
rdistri<- function(n, distri){
  UseMethod("rdistri", distri)
}

#' @rdname numericDistri
#' @export
rdistri.numericDistri<- function(n, distri){
  if(attributes(distri)$logP){
    distri$p<- exp(distri$p)
  }
  sample(distri$x, n, replace=TRUE, prob=distri$p)
}


