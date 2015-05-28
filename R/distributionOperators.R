## Sum ----
distriSum<- function(x, y){
  if (!inherits(x, "numericDistri") | !inherits(y, "numericDistri")) stop("parameter must be numericDistri objects")
  res<- .External("distrisum", x$x, x$p, y$x, y$p, max(x$x) + max(y$x))
  
  res<- data.frame(x=0:(length(res) - 1), p=res)
  infiniteDomain<- inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport")
  attributes(res)$p.omitted<- ifelse(infiniteDomain, 1 - sum(res$p), 0)
  attributes(res)$parameters<- list(x=attributes(x)$parameters, y=attributes(y)$parameters)
  class(res)<- "sumOfDistri"
  if (infiniteDomain) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}
## Test
# gctorture(TRUE) # to detect problems on memory management (very slow)
# distri<- distriBinom(5, .5)
# sapply(1:50, function(x) sum(distriSum(distri, distri)$p))

## TODO: + doesn't work
`+.numericDistri`<- function(x, y){
  return(distriSum(x, y))
}

## R implementation
# distriSum<- function(x, y){
#   names(x)<- c("x", "px")
#   names(y)<- c("y", "py")
#   res<- merge(x, y, all=TRUE)
#   res<- data.frame(x=res$x + res$y, p=res$px * res$py)
#   res<- by(res$p, res$x, sum)
#   res<- data.frame(x=as.numeric(names(res)), p=as.numeric(res))
#   infiniteDomain<- inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport")
#   attributes(res)$p.omitted<- ifelse(infiniteDomain, 1 - sum(res$p), 0)
#   attributes(res)$parameters<- list(x=attributes(x)$parameters, y=attributes(y)$parameters)
#   class(res)<- "sumOfDistri"
#   if (infiniteDomain) class(res)<- c(class(res), "infiniteSuport")
#   class(res)<- c(class(res), "numericDistri", "data.frame")
#   
#   return (res)
# }



## Difference ----
## R implementation
distriDiff<- function(x, y){
  names(x)<- c("x", "px")
  names(y)<- c("y", "py")
  res<- merge(x, y, all=TRUE)
  res<- data.frame(x=res$x - res$y, p=res$px * res$py)
  res<- by(res$p, res$x, sum)
  res<- data.frame(x=as.numeric(names(res)), p=as.numeric(res))
  infiniteDomain<- inherits(x, "infiniteSuport") | inherits(y, "infiniteSuport")
  attributes(res)$p.omitted<- ifelse(infiniteDomain, 1 - sum(res$p), 0)
  attributes(res)$parameters<- list(x=attributes(x)$parameters, y=attributes(y)$parameters)
  class(res)<- "diffOfDistri"
  if (infiniteDomain) class(res)<- c(class(res), "infiniteSuport")
  class(res)<- c(class(res), "numericDistri", "data.frame")
  
  return (res)
}

`-.numericDistri`<- function(x, y){
  return(distriDiff(x, y))
}

