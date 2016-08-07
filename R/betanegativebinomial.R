#' Beta negative binomial probability distribution
#' @name betanbinom
#' @param x 
#' @param q 
#' @param p 
#' @param n 
#' @param size 
#' @param shape1 
#' @param shape2 
#' @param mean 
#' @param variance 
#' @param log 
#'
#' @useDynLib LHR
NULL

#' @describeIn betanbinom
#' @export
dbetanbinom <- function(x, size, shape1, shape2, mean, variance, log = FALSE){
  if (!missing(mean) & !missing(variance)) {
    if (!missing(shape1) | !missing(shape2)) 
      stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
    betaPar<- fbetanbinom(size, mean, variance)
    .External("actuar_do_dpq", "dbetanbinom", x, size, betaPar[[1]], betaPar[[2]], log)
  }
  else .External("actuar_do_dpq", "dbetanbinom", x, size, shape1, shape2, log)
}

#' @describeIn betanbinom
#' @export
pbetanbinom <- function(q, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE){
  if (!missing(mean) & !missing(variance)) {
    if (!missing(shape1) | !missing(shape2)) 
      stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
    betaPar<- fbetanbinom(size, mean, variance)
    .External("actuar_do_dpq", "pbetanbinom", q, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
  }
  else .External("actuar_do_dpq", "pbetanbinom", q, size, shape1, shape2, lower.tail, log.p)
}

#' @describeIn betanbinom
#' @export
qbetanbinom <- function(p, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE){
  if (!missing(mean) & !missing(variance)) {
    if (!missing(shape1) | !missing(shape2))
      stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
    betaPar<- fbetanbinom(size, mean, variance)
    .External("actuar_do_dpq", "qbetanbinom", p, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
  }
  else .External("actuar_do_dpq", "qbetanbinom", p, size, shape1, shape2, lower.tail, log.p)
}

#' @describeIn betanbinom
#' @export
rbetanbinom <- function(n, size, shape1, shape2, mean, variance){
  if (!missing(mean) & !missing(variance)) {
    if (!missing(shape1) | !missing(shape2)) 
      stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
    betaPar<- fbetanbinom(size, mean, variance)
    .External("actuar_do_random", "rbetanbinom", n, size, betaPar[[1]], betaPar[[2]])
  }
  else .External("actuar_do_random", "rbetanbinom", n, size, shape1, shape2)
}

#' @describeIn betanbinom
#' @export
sbetanbinom <- function(size, shape1, shape2){
  mean<- size*shape2/(shape1-1)
  mean[shape1 <= 1]<- Inf
  
  var<- size*(shape1+size-1)*shape2*(shape1+shape2-1)/((shape1-2)*((shape1-1)^2))
  var[shape1 <= 2]<- Inf
  
  if (size < 0 || shape1 < 0 || shape2 < 0){
  	mean[which(size < 0 | shape1 < 0 | shape2 < 0)]<- NaN
  	var[which(is.na(mean))]<- NaN
  	warning("NaNs produced (parameters out of domain)")
  }
  
  return (list(mean=mean, var=var))
}

#' @describeIn betanbinom
#' @export
fbetanbinom <- function(size, mean, var){ # return(data.frame(shape1,shape2))
# Hi ha restriccions en l'espai mean ~ var. alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mean=size*shape2/(shape1-1), var= size*(shape1+size-1)*shape2*(shape1+shape2-1)/((shape1-2)*((shape1-1)^2))], [shape1,shape2]);
  shape1<- (2 * size * var + mean * size^2 + (mean^2 - mean) * size - mean^2) / (size * var - mean * size - mean^2)
  shape2<- (mean * var + mean^2 * size + mean^3) / (size * var - mean * size - mean^2)
  shape1[which(shape1 <= 2 | shape2 <= 0)]<- NaN
  shape2[which(is.na(shape1))]<- NaN
  
  if (size < 0 || mean < 0 || var < 0){
  	shape1[which(size < 0 | mean < 0 | var < 0)]<- NaN
  	shape2[which(is.na(shape1))]<- NaN
  	warning("NaNs produced (parameters out of domain)")
  }

  return (data.frame(shape1, shape2))
}