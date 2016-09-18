##TODO: check ghyper {SuppDists} for documentation pattern
# s* summary functions return the mean and variance of the probability distribution
# f* find functions return the parameters of the probability distribution for the disired mean and var

## Beta distribution
#' Beta moments
#'
#' @param shape1 
#' @param shape2 
#'
#' @return
#' @export
#'
#' @examples
sbeta<- function(shape1, shape2){
  mean<- shape1 / (shape1 + shape2)
  var<- shape1 * shape2 / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))
  if (all(is.na(mean))) return (data.frame(mean=mean, var=mean))
                           
  if (shape1 < 0 || shape2 < 0){
  	mean[which(shape1 < 0 | shape2 < 0)]<- NaN
  	var[which(is.na(mean))]<- NaN
  	warning("NaNs produced: parameters must be > 0")
  }
  
  return (data.frame(mean=mean, var=var))
}

#' Find beta parameters
#'
#' @param mean 
#' @param var 
#'
#' @return
#' @export
#'
#' @examples
fbeta<- function(mean, var="max"){
# Restrictions on the parameter space: max(var) = mean - mean^2 
#   (max var = 0.25 & mean= 0.5)
#   alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mean= shape1/(shape1+shape2) , var= shape1*shape2/((shape1+shape2)^2 * (shape1+shape2+1))], [shape1,shape2]);
  if (var[1] == "max") var<- maxVarBeta(mean) - 0.0000000000001
  shape1<-  -(mean * var + mean^3 - mean^2) / var
  shape2<- ((mean-1) * var + mean^3 - 2 * mean^2 + mean) / var
  
  shape1[which(shape1 <= 0 | shape2 <= 0)]<- NaN
  
  if (mean < 0 || mean > 1){
  	shape1[which(mean < 0 | mean > 1)]<- NaN
  	warning("Mean parameter must be between 0 and 1")
  }
  
  shape2[which(is.na(shape1))]<- NaN
      
  if (any(is.nan(shape1) & var != 0))
  	warning("NaNs produced (parameters out of domain)")

  return (data.frame(shape1, shape2))
}

maxVarBeta<- function(mean) mean - mean^2 # max var is the superior limit (not included)

## Binomial distribution
#' Binomial distribution summary
#'
#' @param size 
#' @param prob 
#'
#' @return
#' @export
#'
#' @examples
sbinom<- function(size, prob){
  mean<- size * prob
  var<- size * prob * (1 - prob)

  if (size < 0 || prob < 0 || prob > 1){
  	mean[which(size < 0 | prob < 0 | prob > 1)]<- NaN
  	var[which(size < 0 | prob < 0 | prob > 1)]<- NaN
  	warning("NaNs produced (parameters out of domain)")
  }
  
  return (data.frame(mean=mean, var=var))
}

## Negative binomial distribution
#' Negative binomial distribution summary
#'
#' @param size 
#' @param prob 
#'
#' @return
#' @export
#'
#' @examples
snbinom<- function(size, prob){
  mean<- prob*size/(1-prob)
  var<- size*prob/((1-prob)^2)

  if (size < 0 || prob < 0 || prob > 1){
  	mean[which(size < 0 | prob < 0 | prob > 1)]<- NaN
  	var[which(size < 0 | prob < 0 | prob > 1)]<- NaN
  	warning("NaNs produced (parameters out of domain)")
  }
  
  return (data.frame(mean=mean, var=var))
}


## Poisson distribution
#' Poisson distribution summary
#'
#' @param lambda 
#'
#' @return
#' @export
#'
#' @examples
spois<- function(lambda){
  mean<- lambda
  var<- lambda
  
  if (any(negL<- lambda < 0)){
    mean[negL]<- NaN
    var[negL]<- NaN
    warning("NaNs produced (parameter out of domain [ lambda > 0 ]).")
  }
  
  return (data.frame(mean=mean, var=var))
}


## Beta negative binomial and Beta binomial
# see betabinomial.R and betanegativebinomial.R