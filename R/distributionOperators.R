## Sum ----

# @export
distriSum<- function(x, y){
  if (!inherits(x, "numericDistri") | !inherits(y, "numericDistri"))
    stop("Parameters must be numericDistri objects.")
    
  if (attributes(x)$logP != attributes(y)$logP){
    y<- logP(y, logP=attributes(x)$logP)
    message("y probabilities transformed to the same scale than x. Use logP(x, log=T/F) to change it.")
  }
  
  log<- attributes(x)$logP
  minRes<- min(x$x) + min(y$x)
  maxRes<- max(x$x) + max(y$x)
  vals<- minRes:maxRes
  
  res<- .External("distrisum", x$x, x$p, y$x, y$p, log, minRes, maxRes, length(vals))
  res<- data.frame(x=vals, p=res)
  
  infiniteDomain<- inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport")
  if (log){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(x=attributes(x)$parameters, y=attributes(y)$parameters)
  attributes(res)$logP<- log
  class(res)<- "sumOfDistri"
  if (infiniteDomain) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}
## Test
# gctorture(TRUE) # to detect problems on memory management (very slow)
# distri<- distriBinom(5, .5)
# sapply(1:50, function(x) sum(distriSum(distri, distri)$p))

#' @rdname numericDistri
#' @export
"+.numericDistri"<- function(x, y){
  return(distriSum(x, y))
}


## Difference ----

## R implementation
# @export
distriDiff<- function(x, y){
  if (!inherits(x, "numericDistri") | !inherits(y, "numericDistri")) stop("Parameters must be numericDistri objects.")
  if (attributes(x)$logP | attributes(y)$logP) stop("Parameters must have probabilities on natural scalse. Use logP(x, log=F) to change it. Perhaps I will implement it some day.")
  log<- attributes(x)$logP

  names(x)<- c("x", "px")
  names(y)<- c("y", "py")
  res<- merge(x[x$px > 0, ], y[y$py > 0, ], all=TRUE) ## Remove values with p = 0
  res<- data.frame(x=res$x - res$y, p=res$px * res$py)
  res<- by(res$p, res$x, sum)
  res<- data.frame(x=as.numeric(names(res)), p=as.numeric(res))

  infiniteDomain<- inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport")
  if (log){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(x=attributes(x)$parameters, y=attributes(y)$parameters)
  attributes(res)$logP<- log
  class(res)<- "diffOfDistri"
  if (infiniteDomain) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}

#' @rdname numericDistri
#' @export
"-.numericDistri"<- function(x, y){
  return(distriDiff(x, y))
}

## Scalar product ----

# @export
distriScalarProd<- function(distri, x){
  if (x == 1) return (distri)
  
  distri$x<- distri$x * x
  
#   domain<- 0:max(distri$x)
#   domain<- domain[!domain %in% distri$x]
#   domain<- data.frame(x=domain, p=0)
#   
#   distri<- rbind(distri, domain)
#   distri<- distri[order(distri$x),]
#   rownames(distri)<- NULL
  
  attributes(distri)$parameters<- c(list(scalarProd=x), attributes(distri)$parameters)
  class(distri)<- c("scalarProdDistri", class(distri))
  
  return(distri)
}

#' @rdname numericDistri
#' @param x a positive integer
#' @export
"*.numericDistri"<- function(distri, x){
  return(distriScalarProd(distri, x))
}
