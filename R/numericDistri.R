# Functions for compound discrete probability distributions resolved numerically
# numericDistri class

#' Numeric representation of discrete distributions
#' 
#' @name numericDistri
NULL

## TODO: call distriBeta* or distri* according to p parameter (p = c(shape1, shape2) | p)
## TODO: optimize for prob=1 & prob=0

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
  if (size < 0)
    stop("Negative size parameter not allowed in Binomial distributions.")
  
  res<- data.frame(x=0:size, p=dbinom(x=0:size, size, prob, log=logP))
  attributes(res)$p.omitted<- 0
  attributes(res)$parameters<- list(size=size, prob=prob)
  attributes(res)$support<- c(0, size)
  attributes(res)$logP<- logP
  class(res)<- c("binom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# size: numericDistri object which compound a binomial distribution as a size parameter.
# p: probability
#' @rdname numericDistri
#' @export
distriBinom.numericDistri<- function(size, prob, logP=FALSE){
  if (min(size$x) < 0)
    stop("Negative size parameter not allowed in Binomial distributions.")
  
  if (logP != attributes(size)$logP){ # Transform p to log scale if necessary
    size<- logP(size, logP=logP)
  }
  
  if (logP){
    size<- size[size$p > log(getOption("numericDistri.minP")),, drop=FALSE]
  }else{
    size<- size[size$p > getOption("numericDistri.minP"),, drop=FALSE]
  }
  
  if (nrow(size) == 0)
    stop("No value with probability > getOption(\"numericDistri.minP\").\n",
         "\tUse: options(list(numericDistri.minP=1e-100) to change it.\n")
  
  maxRes<- max(size$x)

  res<- .External("binomialCompound", size$x, size$p, prob, logP, maxRes)
  res<- data.frame(x=0:maxRes, p=res)
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(size=class(size)[1], prob=prob)
  attributes(res)$support<- c(0, max(attributes(size)$support))
  attributes(res)$logP<- logP
  
  class(res)<- "compoundBinom"
  
  if (inherits(size, "infiniteSuport"))
    class(res)<- c(class(res), "infiniteSuport")
    
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
  if (size < 0)
    stop("Negative size parameter not allowed in Beta Binomial distributions.")
  
  res<- data.frame(x=0:size, p=dbetabinom(x=0:size, size, shape1, shape2, log=logP))
  attributes(res)$p.omitted<- 0
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var)
  attributes(res)$support<- c(0, size)
  attributes(res)$logP<- logP
  class(res)<- c("betaBinom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}


# size: numericDistribution object which compound a binomial distribution as a size parameter.
# p: probability
#' @rdname numericDistri
#' @export
distriBetaBinom.numericDistri<- function(size, shape1, shape2, logP=FALSE){
  if (min(size$x) < 0)
    stop("Negative size parameter not allowed in Beta Binomial distributions.")
  
  if (logP != attributes(size)$logP){ # Transform p to log scale if necessary
    size<- logP(size, logP=logP)
  }
  
  if (logP){
    size<- size[size$p > log(getOption("numericDistri.minP")),, drop=FALSE]
  }else{
    size<- size[size$p > getOption("numericDistri.minP"),, drop=FALSE]
  }
  
  if (nrow(size) == 0)
    stop("No value with probability > getOption(\"numericDistri.minP\").\n",
         "\tUse: options(list(numericDistri.minP=1e-100) to change it.\n")
  
  maxRes<- max(size$x)
  
  res<- .External("BetaBinomialCompound", size$x, size$p, shape1, shape2, logP, maxRes)
  res<- data.frame(x=0:maxRes, p=res)
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=class(size)[1], prob=prob$mean, var=prob$var)
  attributes(res)$support<- c(0, max(attributes(size)$support))
  attributes(res)$logP<- logP
  
  class(res)<- "compoundBetaBinom"
  
  if (inherits(size, "infiniteSuport"))
    class(res)<- c(class(res), "infiniteSuport")
    
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
distriNegBinom<- function(size, prob, mu, logP=FALSE, minP, maxPomitted=0.01, maxX){
  UseMethod("distriNegBinom")
}

#' @rdname numericDistri
#' @export  
distriNegBinom.numeric<- function(size, prob, mu, logP=FALSE, minP, maxPomitted=0.01, maxX){
  if (missing(prob) & !missing(mu)) prob<- size / (size+mu) #alternative parametrization ?dnbinom
  if (missing(maxX)){
    maxX<- qnbinom(maxPomitted, size=size, prob=prob, lower.tail=FALSE, log.p=logP)
  }
  res<- data.frame(x=0:maxX, p=dnbinom(x=0:maxX, size=size, prob=prob, log=logP))
  
  if (!missing(minP)) res<-  res[res$p > minP,]
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(size=size, prob=prob)
  attributes(res)$support<- c(0, Inf)
  attributes(res)$logP<- logP
  class(res)<- c("nbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}

## TODO: distriNegBinom.numericDistri ----

## Beta negative binomial distribution ----
#' @rdname numericDistri
#' @export
distriBetaNegBinom<- function(size, shape1, shape2, logP=FALSE, minP, maxPomitted=0.01, maxX){
  if (missing(maxX)){
    maxX<- qbetanbinom(maxPomitted, size=size, shape1=shape1, shape2=shape2, lower.tail=FALSE, log.p=logP)
  }
  res<- data.frame(x=0:maxX, p=dbetanbinom(x=0:maxX, size=size, shape1=shape1, shape2=shape2, log=logP))
  
  if (!missing(minP)) res<-  res[res$p > minP,]
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var)
  attributes(res)$support<- c(0, Inf)
  attributes(res)$logP<- logP
  class(res)<- c("betaNbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}

#TODO distriBetaNegBinom.numericDistri ----


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
distriPois<- function(lambda, logP=FALSE, minP, maxPomitted=0.0001){
  UseMethod("distriPois")
}

#' @rdname numericDistri
#' @export
distriPois.numeric<- function(lambda, logP=FALSE, minP, maxPomitted=0.0001){
  domain<- 0:qpois(1 - maxPomitted, lambda=lambda)
  res<- data.frame(x=domain, p=dpois(x=domain, lambda, log=logP))
  
  if (!missing(minP)) res<-  res[res$p > minP,]
  
  if (logP){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  
  attributes(res)$parameters<- list(lambda=lambda)
  attributes(res)$support<- c(0, Inf)
  attributes(res)$logP<- logP
  class(res)<- c("pois", "infiniteSuport","numericDistri", "data.frame")
  return (res)
}


## TODO: distriPois.numericDistri (src/numericDistribution.c) ----
# @rdname numericDistri
# @export
distriPois.numericDistri<- function(lambda, logP=FALSE, minP, maxPomitted=0.0001){
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
  attributes(res)$support<- c(0, Inf)
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
print.numericDistri<- function(x, maxRows=20, ...){
  cat("\t", class(x)[1], "object\n")
  cat("Theoretical support: [", paste(attributes(x)$support, collapse=", "), "]\n", sep="")
  cat("Probability omitted:", attributes(x)$p.omitted, "\n")
  if (attributes(x)$logP) cat("p in Log probability scale\n")
  
  cat("Parameters:\n")
  utils::str(attributes(x)$parameters)
  
  cat("\n")
  print(utils::head(data.frame(x), n=maxRows), ...)
  if (nrow(x) > maxRows) cat("\t ...\t", nrow(x) - maxRows, "rows omited.\n")
  
  invisible(x)
}

#' @rdname numericDistri
#' @export
summary.numericDistri<- function(object, ...){
  res<- cbind(sdistri(object), t(quantile(object)))

  at<- attributes(object)
  attributes(res)$p.omitted<- at$p.omitted
  attributes(res)$parameters<- at$parameters
  attributes(res)$support<- at$support
  attributes(res)$logP<- at$logP
  
  class(res)<- c("summary.numericDistri", class(res))

  return (res)
}

#' @export
#' @S3method print summary.numericDistri
print.summary.numericDistri<- function(x, ...){
  cat("\t", class(x)[1], "object\n")
  cat("Theoretical support: [", paste(attributes(x)$support, collapse=", "), "]\n", sep="")
  cat("Probability omitted:", attributes(x)$p.omitted, "\n")
  if (attributes(x)$logP) cat("p in Log probability scale\n")
  
  cat("Parameters:\n")
  utils::str(attributes(x)$parameters)

  cat("\n")
  print(as.data.frame(x), row.names=FALSE, ...)
}

#' @rdname numericDistri
#' @export
plot.numericDistri<- function(x, cum=FALSE, type=ifelse(cum, "s", "p"), ...){
  if (cum){
    x<- cumP(x)
    x$p<- x$cump
  }

  graphics::plot.default(x, type=type, ...)
}

#' @rdname numericDistri
#' @export
points.numericDistri<- function(x, cum=FALSE, type=ifelse(cum, "s", "p"), ...){
  if (cum){
    x<- cumP(x)
    x$p<- x$cump
  }
  
  graphics::points.default(x, type=type, ...)
}

## Stats ----
#' @rdname numericDistri
#' @export
mean.numericDistri<- function(x, ...){
  distri<- logP(x, logP=FALSE)
  distri<- distri[distri$p > 0, ]
  
  res<- stats::weighted.mean(distri$x, distri$p, ...)

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
  distri<- distri[distri$p > 0, ]
  
  sum(distri$p * (distri$x - mean(distri))^2)
}

#' @rdname numericDistri
#' @export
median.numericDistri<- function(x, na.rm=FALSE, ...){
  qdistri(.5, x)
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
  if ("cump" %in% names(x))
    return(x)
  
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
quantile.numericDistri<- function(x, na.rm=FALSE, probs=seq(0, 1, 0.25), ...){
  res<- qdistri(probs, x)
  names(res)<- paste0(round(probs * 100, digits=1), "%")
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
    
    # value of the random variable such that the probability of the variable being less than or equal to that value
    sel<- which.min(abs(diffs))
    
    if (diffs[sel] > 0 & sel > 1){
      sel<- sel - 1
    }
    sel
  })
  
  res<- distri$x[sel]
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
