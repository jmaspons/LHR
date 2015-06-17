# Functions for compound discrete probability distributions resolved numerically
# numericDistri class 
setOldClass("numericDistri")

## TODO: call distriBeta* or distri* according to p parameter (p = c(shape1, shape2) | p)

## Binomial distribution ----
distriBinom<- function(...){
  UseMethod("distriBinom")
}

distriBinom.numeric<- function(size, prob, log=FALSE){
  res<- data.frame(x=0:size, p=dbinom(x=0:size, size, prob, log))
  attributes(res)$p.omitted<- 0
  attributes(res)$parameters<- list(size=size, prob=prob, log=log)
  class(res)<- c("binom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# size: numericDistri object which compound a binomial distribution as a size parameter.
# p: probability
distriBinom.numericDistri<- function(size, prob, log=FALSE){
  res<- .External("binomialCompound", size$p, prob, log)
  res<- data.frame(x=0:(length(res) - 1), p=res)
  
  attributes(res)$p.omitted<- 1 - sum(res$p)
  attributes(res)$parameters<- list(size=c(list(distribution=class(size)[1]), attributes(size)$parameters), prob=prob, log=log)
  class(res)<- "compoundBinom"
  if (inherits(size, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}

## Beta binomial distribution ----
distriBetaBinom<- function(...){
  UseMethod("distriBetaBinom")
}

distriBetaBinom.numeric<- function(size, shape1, shape2, log=FALSE){
  res<- data.frame(x=0:size, p=dbetabinom(x=0:size, size, shape1, shape2, log=log))
  attributes(res)$p.omitted<- 0
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var, log=log)
  class(res)<- c("betaBinom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# size: numericDistribution object which compound a binomial distribution as a size parameter.
# p: probability
distriBetaBinom.numericDistri<- function(size, shape1, shape2, log=FALSE){
  res<- .External("BetaBinomialCompound", size$p, shape1, shape2, log)
  res<- data.frame(x=0:(length(res) - 1), p=res)
  
  attributes(res)$p.omitted<- 1 - sum(res$p)
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=c(list(distribution=class(size)[1]), attributes(size)$parameters), prob=prob$mean, var=prob$var, log=log)
  class(res)<- "compoundBetaBinom"
  if (inherits(size, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}

## Negative binomial distribution ----
distriNegBinom<- function(...){
  UseMethod("distriNegBinom")
}
  
distriNegBinom.numeric<- function(size, prob, mu, log=FALSE, maxPomitted=0.01, maxX){
  if (missing(prob) & !missing(mu)) prob<- size / (size+mu) #alternative parametrization ?dnbinom
  if (missing(maxX)){
    maxX<- qnbinom(maxPomitted, size=size, prob=prob, lower.tail=FALSE, log.p=log)
  }
  res<- data.frame(x=0:maxX, p=dnbinom(x=0:maxX, size=size, prob=prob, log=log))
  attributes(res)$p.omitted<- 1 - sum(res$p)
  attributes(res)$parameters<- list(size=size, prob=prob, log=log)
  class(res)<- c("nbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}

#TODO distriNegBinom.numericDistri

## Beta negative binomial distribution ----
distriBetaNegBinom<- function(size, shape1, shape2, log=FALSE, maxPomitted=0.01, maxX){
  if (missing(maxX)){
    maxX<- qbetanbinom(maxPomitted, size=size, shape1=shape1, shape2=shape2, lower.tail=FALSE, log.p=log)
  }
  res<- data.frame(x=0:maxX, p=dbetanbinom(x=0:maxX, size=size, shape1=shape1, shape2=shape2, log=log))
  attributes(res)$p.omitted<- 1 - sum(res$p)
  prob<- sbeta(shape1=shape1, shape2=shape2)
  attributes(res)$parameters<- list(size=size, prob=prob$mean, var=prob$var, log=log)
  class(res)<- c("betaNbinom", "infiniteSuport", "numericDistri", "data.frame")
  return (res)
}


## Generic methods for numericDistri class ----
print.numericDistri<- function(distri){
  cat("\t", class(distri)[1], "distribution\n")
  if (attributes(distri)$p.omitted > 0) cat("Probability omitted: ", attributes(distri)$p.omitted, "\n")
  cat("Parameters:\n")
  str(attributes(distri)$parameters)
  cat("\n")
  print(head(data.frame(distri), n=15))
  if (nrow(distri) > 15) cat("\t ...\t", nrow(distri) - 15, "rows omited.\n")
}

summary.numericDistri<- function(distri){
  cat("\t", class(distri)[1], "distribution\n")
  if (attributes(distri)$p.omitted > 0) cat("Probability omitted: ", attributes(distri)$p.omitted, "\nParameters:\n")
  str(attributes(distri)$parameters)
  cat("\n")
  res<- sdistri(distri)
  print(res)
  invisible(res)
}


plot.numericDistri<- function(distri, cum=FALSE, ...){
  if (cum){
    distri$p<- cumsum(distri$p)
  }
  type<- ifelse(cum, "s", "p")
  plot.default(distri, type=type, ...)
}


## Stats ----
mean.numericDistri<- function(distri){
  weighted.mean(distri$x, distri$p)
}

var.default<- var

var<- function(distri, ...){
  UseMethod("var")
}

var.numericDistri<- function(distri){
  sum(distri$p * (distri$x - mean(distri))^2)
}

sdistri<- function(distri){
  UseMethod("sdistri")
}

sdistri.numericDistri<- function(distri){
  meanD<- mean(distri)
  varD<- sum(distri$p * (distri$x - meanD)^2) #var(distri)
  GmeanD<- G(meanD, varD)
  
  return (data.frame(Gmean=GmeanD, mean=meanD, var=varD))
}

cumsum.numericDistri<- function(distri){
  distri$cump<- cumsum(distri$p)
  return(distri)
}


# random deviates from a numericDistribution
rdistri<- function(n, distri){
  UseMethod("rdistri", distri)
}

rdistri.numericDistri<- function(n, distri){
  sample(distri$x, n, replace=TRUE, prob=distri$p)
}



