## Sum ----

# @export
distriSum<- function(x, y){
  if (!inherits(x, "numericDistri") | !inherits(y, "numericDistri")){
    if (is.numeric(x)){
      y$x<- y$x + x
      attributes(y)$support<- attributes(y)$support + x
      class(y)<- c("numericDistri", "sumOfDistri")
      if (inherits(y, "infiniteSuport"))
        class(y)<- c(class(y), "infiniteSuport")
      
      class(y)<- c(class(y), "data.frame")
      return(y)
    }
    if (is.numeric(y)){
      x$x<- x$x + y
      attributes(x)$support<- attributes(x)$support + y
      class(x)<- c("numericDistri", "sumOfDistri")
      if (inherits(x, "infiniteSuport"))
        class(x)<- c(class(x), "infiniteSuport")
      
      class(x)<- c(class(x), "data.frame")
      return(x)
    }
    stop("Parameters must be numericDistri objects or numeric. ")
  }
  
  if (attributes(x)$logP != attributes(y)$logP){
    y<- logP(y, logP=attributes(x)$logP)
    warning("y probabilities transformed to the same scale than x. Use logP(x, log=TRUE/FALSE) to change it.")
  }
  
  if (attributes(x)$logP){
    x<- x[x$p > -Inf,, drop=FALSE]
    y<- y[y$p > -Inf,, drop=FALSE]
  }else{
    x<- x[x$p > 0,, drop=FALSE]
    y<- y[y$p > 0,, drop=FALSE]
  }
  
  log<- attributes(x)$logP
  minRes<- min(x$x) + min(y$x)
  maxRes<- max(x$x) + max(y$x)
  vals<- minRes:maxRes
  
  res<- .External("distrisum", x$x, x$p, y$x, y$p, log, minRes, length(vals))
  res<- data.frame(x=vals, p=res)
  
  if (log){
    attributes(res)$p.omitted<- 1 - sum(exp(res$p))
  }else{
    attributes(res)$p.omitted<- 1 - sum(res$p)
  }
  attributes(res)$parameters<- list(x=class(x)[1], y=class(y)[1])
  attributes(res)$support<- attributes(x)$support + attributes(y)$support
  attributes(res)$logP<- log
  
  class(res)<- c("sumOfDistri", "numericDistri")
  if (inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport"))
    class(res)<- c(class(res), "infiniteSuport")
    
  class(res)<- c(class(res), "data.frame")
  
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
  if (inherits(y, "numericDistri")){
    y$x<- -y$x
    attributes(y)$support<- range(- attributes(y)$support)
  }else{
    y<- -y
  }
  
  distriSum(x, y)
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
  
  attributes(distri)$parameters<- c(list(scalarProd=x), attributes(distri)$parameters)
  attributes(distri)$support<- attributes(distri)$support * x
  class(distri)<- c("scalarProdDistri", class(distri))
  
  return(distri)
}

#' @rdname numericDistri
#' @param x a positive integer (not checked).
#' @export
"*.numericDistri"<- function(distri, x){
  return(distriScalarProd(distri, x))
}


## Add negative values to 0 ----
## TODO: logP
#' Title
#'
#' @param distri 
#'
#' Useful to simulate natural values (e.g. populations)
#' @return
#' @export
#'
#' @examples
neg2zeroP<- function(distri){
  if (attributes(distri)$logP){
    distri<- logP(distri, logP=FALSE)
    warning("Probabilities transformed to non Log scale.")
  }
  
  if (!any(distri$x == 0)){
    distri<- rbind(distri, data.frame(x=0, p=0))
    distri<- distri[order(distri$x), ]
  }
  
  negP<- distri$p[distri$x < 0]
  distri<- distri[distri$x >= 0, ]
  
  distri$p[distri$x == 0]<- distri$p[distri$x == 0] + sum(negP)
  maxSupport<- ifelse(max(attributes(distri)$support) > 0, max(attributes(distri)$support), 0)
  attributes(distri)$support<- c(0, maxSupport)
  
  class(distri)<- c("numericDistri", "data.frame")
  
  return(distri)
}
