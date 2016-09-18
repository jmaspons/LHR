# Functions for compound discrete probability distributions resolved numerically
# numericDistri class

#' Numeric representation of discrete distributions
#' 
#' @name numericDistri
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


## Poisson distribution ----
#' @rdname numericDistri
#'
#' @inheritParams distriBinom
#' @param lambda 
#' @param logP 
#' @param minP probabilities smaller than \code{minP} are discarded. Useful to reduce memory usage.
#'
#' @return
#' @importFrom stats dpois dpois qpois qpois
#' @export
distriPois<- function(lambda, logP=FALSE, minP=.Machine$double.eps){
  UseMethod("distriPois")
}

#' @rdname numericDistri
#' @export
distriPois.numeric<- function(lambda, logP=FALSE, minP=.Machine$double.eps){
  domain<- max(50, lambda * 3)
  domain<- 0:domain
  res<- data.frame(x=domain, p=dpois(x=domain, lambda, log=logP))
  
  if (!missing(minP)) res<-  res[res$p > minP,]
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  
  attributes(res)$parameters<- list(lambda=lambda)
  attributes(res)$logP<- logP
  class(res)<- c("pois", "infiniteSuport","numericDistri", "data.frame")
  return (res)
}


## TODO
# @rdname numericDistri
# @export
distriPois.numericDistri<- function(lambda, logP=FALSE, minP=.Machine$double.eps){
  stop("Compound Poisson distribution not implemented.")
  ## TODO
  
  if (logP != attributes(lambda)$logP){ # Transform p to log scale if necessary
    lambda<- logP(lambda, logP=logP)
  }
  
  res<- .External("poissonCompound", lambda$x, lambda$p, lambda, logP, max(lambda$x))
  res<- data.frame(x=0:(length(res) - 1), p=res)
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(lambda=c(list(distribution=class(lambda)[1]), attributes(lambda)$parameters))
  attributes(res)$logP<- logP
  
  class(res)<- c("compoundPois", "infiniteSuport", "numericDistri", "data.frame")
  
  return (res)
}


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

# @export
# print<- function(x, ...) base::print(x, ...)

#' @export
#' @S3method print numericDistri
print.numericDistri<- function(x, ...){
  cat("\t", class(x)[1], "distribution\n")
  cat("Probability omitted: ", attributes(x)$p.omitted, "\n")
  if (attributes(x)$logP) cat("p in Log probability scale\n")
  cat("Parameters:\n")
  utils::str(attributes(x)$parameters)
  cat("\n")
  print(utils::head(data.frame(x), n=15), ...)
  if (nrow(x) > 15) cat("\t ...\t", nrow(x) - 15, "rows omited.\n")
  
  invisible(x)
}

#' @rdname numericDistri
#' @export
summary.numericDistri<- function(object, ...){
  res<- cbind(sdistri(object), t(quantile(object)))

  at<- attributes(object)
  attributes(res)$p.omitted<- at$p.omitted
  attributes(res)$parameters<- at$parameters
  attributes(res)$logP<- at$logP
  
  class(res)<- c("summary.numericDistri", class(res))

  return (res)
}

#' @export
#' @S3method print summary.numericDistri
print.summary.numericDistri<- function(x, ...){
  cat("Probability omitted: ", attributes(x)$p.omitted, "\nParameters:\n")
  
  if (attributes(x)$logP) cat("p in Log probability scale\n")
  
  utils::str(attributes(x)$parameters)
  
  cat("\n")
  print(as.data.frame(x), row.names=FALSE, ...)
}

#' @rdname numericDistri
#' @export
plot.numericDistri<- function(x, cum=FALSE, ...){
  if (cum){
    x<- cumP(x)
    x$p<- x$cump
  }
  type<- ifelse(cum, "s", "p")
  graphics::plot.default(x, type=type, ...)
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


# @exportMethod var
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
sdistri<- function(x){
  UseMethod("sdistri")
}

#' @rdname numericDistri
#' @export
sdistri.numericDistri<- function(x){
  distri<- logP(x, logP=FALSE)

  meanD<- mean(distri)
  varD<- sum(distri$p * (distri$x - meanD)^2) #var(distri)
  GmeanD<- G(meanD, varD)
  
  return (data.frame(Gmean=GmeanD, mean=meanD, var=varD))
}


# @exportMethod cumP
#' @export
cumP<- function(x, ...){
  UseMethod("cumP")
}


# @method cumP numericDistri
#' @rdname numericDistri
#' @export
cumP.numericDistri<- function(x, ...){
  if(attributes(x)$logP){
    p<- exp(x$p)
    x$cump<- cumsum(p)
  }else{
    x$cump<- cumsum(x$p)
  }
  
  return(x)
}


#' @rdname numericDistri
#' @export
#' @importFrom stats quantile
quantile.numericDistri<- function(x, na.rm=FALSE, ...){
  x<- cumP(x)
  p<- stats::quantile(x$cump, na.rm=na.rm, ...)
  res<- qdistri(p, x, ...)
  names(res)<- names(p)
  res
}


## Probability distribution ----

#' @rdname numericDistri
#' @export
ddistri<- function(x, distri){
  distri$p[match(x, distri$x)]
}

#' @rdname numericDistri
#' @export
pdistri<- function(q, distri){
  distri<- cumP(distri)
  sel<- sapply(q, function(x){
    diffs<- distri$x - x
    sel<- which.min(abs(diffs))
    
    if (diffs[sel] > 0){
      sel<- sel - 1
    }
    sel
  })
  res<- distri$cump[sel]
  res[sel == 0]<- 0
  res
}

#' @rdname numericDistri
#' @export
qdistri<- function(p, distri){
  distri<- cumP(distri)
  
  if (anyNA(distri$cump)) return(rep(NA_real_, length(p)))
  
  sel<- sapply(p, function(x){
    diffs<- distri$cump - x

    # sel<- which.min(abs(diffs))
    
    sel<- which(abs(diffs) == min(abs(diffs)))
    # Select the x value in the midle if there are more than one value with the minimum distance to p
    if (length(sel) > 1)
      sel<- sel[length(sel) %/% 2]
    
    if (diffs[sel] > 0){
      sel<- sel - 1
    }
    sel
  })
  res<- distri$x[sel]
  res[sel == 0]<- 0
  res[p >= 1]<- max(distri$x)
  res
}

#' Random deviates from a numericDistribution
#' 
#' @export
rdistri<- function(n, x){
  UseMethod("rdistri", x)
}

#' @rdname numericDistri
#' @export
rdistri.numericDistri<- function(n, x){
  if(attributes(x)$logP){
    x$p<- exp(x$p)
  }
  sample(x$x, n, replace=TRUE, prob=x$p)
}