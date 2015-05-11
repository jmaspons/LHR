# Functions for compound discrete probability distributions resolved numerically
# numericDistribution class 
setOldClass("numericDistri")

## TODO: call distriBeta* or distri* according to p parameter (p = c(shape1, shape2) | p)
## TODO: port critical operations to c: distri*.numericDistri

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

# distri: numericDistribution object which compound a binomial distribution as a size parameter.
# p: probability
distriBinom.numericDistri<- function(distri, prob, log=FALSE){
  if (!inherits(distri, "numericDistri")) stop("distri should be a numericDistri object")
  
  sizeDistri<- lapply(1:nrow(distri), function(i){
    x<- distri$x[i]
    p<- distri$p[i]
    res<- data.frame(x=0:x, p=dbinom(0:x, size=x, prob=prob, log=log) * p)
    #     class(res)<- c("numericDistri")
    return (res)
  })
  sizeDistri<- do.call("rbind", sizeDistri)
  res<- by(sizeDistri$p, sizeDistri$x, sum)
  res<- data.frame(x=as.numeric(names(res)), p=as.numeric(res))
  attributes(res)$p.omitted<- 1 - sum(res$p)
  attributes(res)$parameters<- list(size=c(list(distribution=class(distri)[1]), attributes(distri)$parameters), prob=prob, log=log)
  class(res)<- "compoundBinom"
  if (inherits(distri, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
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
  attributes(res)$parameters<- list(size=size, shape1=shape1, shape2=shape2, log=log)
  class(res)<- c("betaBinom", "finiteSuport","numericDistri", "data.frame")
  return (res)
}

# distri: numericDistribution object which compound a binomial distribution as a size parameter.
# p: probability
distriBetaBinom.numericDistri<- function(distri, shape1, shape2){
  if (!inherits(distri, "numericDistri")) stop("distri should be a numericDistri object")
  
  sizeDistri<- lapply(1:nrow(distri), function(i){
    x<- distri$x[i]
    p<- distri$p[i]
    res<- data.frame(x=0:x, p=dbetabinom(0:x, size=x, shape1, shape2) * p)
    #     class(res)<- c("numericDistri")
    return (res)
  })
  sizeDistri<- do.call("rbind", sizeDistri)
  res<- by(sizeDistri$p, sizeDistri$x, sum)
  res<- data.frame(x=as.numeric(names(res)), p=as.numeric(res))
  attributes(res)$p.omitted<- 1 - sum(res$p)
  attributes(res)$parameters<- list(size=c(list(distribution=class(distri)[1]), attributes(distri)$parameters), prob=prob, log=log)
  class(res)<- "compoundBetaBinom"
  if (inherits(distri, "infiniteSuport")) class(res)<- c(class(res), "infiniteSuport")
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
  attributes(res)$parameters<- list(size=size, shape1=shape1, shape2=shape2, log=log)
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
  print(sdistri(distri))
  
}


# random deviates from a numericDistribution
rdistri<- function(n, distri){
  UseMethod("rdistri", distri)
}

rdistri.numericDistri<- function(n, distri){
  sample(distri$x, n, replace=TRUE, prob=distri$p)
}

# stats
sdistri<- function(distri){
  UseMethod("sdistri")
}

sdistri.numericDistri<- function(distri){
  mean<- weighted.mean(distri$x, distri$p)
  var<- sum(distri$p * (distri$x - mean)^2)
  Gmean<- G(mean, var)
  
  return (data.frame(mean, var, Gmean))
}



