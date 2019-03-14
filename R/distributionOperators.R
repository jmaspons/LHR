## Sum ----

# @export
distriSum<- function(x, y){
  if (!inherits(x, "numericDistri") | !inherits(y, "numericDistri")){
    if (is.numeric(x)){
      y$x<- y$x + x
      return(y)
    }
    if (is.numeric(y)){
      x$x<- x$x + y
      return(x)
    }
    stop("At least one parameter must be a numericDistri object.")
  }
  
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
  y$x<- -y$x
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
  
  return(distri)
}
