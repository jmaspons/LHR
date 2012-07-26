
dbetabinom <- function(x, size, shape1, shape2, mean, variance, log = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- fbetabinom(size, mean, variance)
        .External("actuar_do_dpq", "dbetabinom", x, size, betaPar[[1]], betaPar[[2]], log)
    }
    else .External("actuar_do_dpq", "dbetabinom", x, size, shape1, shape2, log)
}

pbetabinom <- function(q, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- fbetabinom(size, mean, variance)
        .External("actuar_do_dpq", "pbetabinom", q, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "pbetabinom", q, size, shape1, shape2, lower.tail, log.p)
}

qbetabinom <- function(p, size, shape1, shape2, mean, variance, lower.tail = TRUE, log.p = FALSE)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- fbetabinom(size, mean, variance)
        .External("actuar_do_dpq", "qbetabinom", p, size, betaPar[[1]], betaPar[[2]], lower.tail, log)
    }
    else .External("actuar_do_dpq", "qbetabinom", p, size, shape1, shape2, lower.tail, log.p)
}

rbetabinom <- function(n, size, shape1, shape2, mean, variance)
{
    if (!missing(mean) & !missing(variance)) {
        if (!missing(shape1) | !missing(shape2)) 
            stop("Beta parameters 'shape1' and 'shape2', and 'mean' and 'variance' both specified")
        betaPar<- fbetabinom(size, mean, variance)
        .External("actuar_do_random", "rbetabinom", n, size, betaPar[[1]], betaPar[[2]])
    }
    else .External("actuar_do_random", "rbetabinom", n, size, shape1, shape2)
}

sbetabinom <- function(size, shape1, shape2)
{
    mean<- size * shape1 / (shape1 + shape2)
    var<- size*shape1*shape2 *(shape1+shape2+size) / (((shape1+shape2)^2) * (shape1+shape2+1))
    
    if (size < 0 || shape1 < 0 || shape2 < 0){
	mean[which(size < 0 | shape1 < 0 | shape2 < 0)]<- NaN
	var[which(is.na(mean))]<- NaN
	warning("NaNs produced (parameters out of domain)")
    }
    
    return (list(mean=mean, var=var))
}


fbetabinom<- function(size, mean, var)
{
# Hi ha restriccions en l'espai mean ~ var: alpha > 1 & beta > 1 -> unimodal
# Maxima: solve([mean=size*shape1/(shape1+shape2) , var= size*shape1*shape2*(shape1+shape2+size)/(((shape1+shape2)^2)*(shape1+shape2+1))], [shape1,shape2]);
    shape1<- -(mean * var - mean^2 * size + mean^3) / (size * var - mean * size + mean^2)
    shape2<- -((size - mean) * var - mean * size^2 + 2 * mean^2 * size - mean^3)/ (size * var - mean * size + mean^2)
    
    shape1[which(shape1 <= 0 | shape2 <= 0)]<- NaN
    shape2[which(is.na(shape1))]<- NaN
  
    if (size < 0 || mean < 0 || var < 0){
	shape1[which(size < 0 | mean < 0 | var < 0)]<- NaN
	shape2[which(is.na(shape1))]<- NaN
	warning("NaNs produced (parameters out of domain)")
    }
  
  return (data.frame(shape1, shape2))
}