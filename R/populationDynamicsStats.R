# pop: matrix(nrow=replicates, ncol=tf+1, dimnames=list(replicate=NULL, t=0:tf)) from models-discreteTime.R
# pop<- mSurvBV.t(broods=2, b=2, breedFail=.5, j=.5, a=.7, N0=20, replicates=100, tf=10)
#' @rdname discretePopSim
#'
#' @param object 
#' @param dt 
#'
#' @return
#' @export
#'
#' @examples
summary.discretePopSim<- function(object, dt=1){
  R<- as.numeric(r(object, dt=dt)) # intrinsic growth rate
  L<- as.numeric(lambda(object, dt=dt)) # lambda
  meanR<- mean(R, na.rm=TRUE)
  varR<- var(R, na.rm=TRUE)
  meanL<- mean(L, na.rm=TRUE)
  varL<- var(L, na.rm=TRUE)
  GR<- G(meanR, varR)
  GL<- G(meanL, varL)
  trends<- trendsProp(object)
  res<- data.frame(trends, GR, meanR, varR, GL, meanL, varL)
  return(res)
}

#' @export
r<- function(...){
  UseMethod("r")  
}

#' @rdname discretePopSim
#' @export
r.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  dN<- t(diff(t(pop))) # dN =  N_t+1 - N_t
  N0<- pop[,-ncol(pop)]
  r<- (dN / dt) / N0
#   names(r)<- colnames(dN) # otherwise it takes the names from N0 (0:(tf-1) instead of 1:tf as does lambda function
  return (r) # intrinsic grow rate (r = dN / dt / N)
}

# r.numericDistri on models-compoundDistributions.R

#' @rdname discretePopSim
#' @export
lambda.discretePopSim<- function(x, dt=1){
  sampleT<- seq(1, ncol(x), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  x<- x[,sampleT]
  
  return (x[,-1] / x[,-ncol(x)]) # lambda = Nt+1 / Nt
}

# lambda.leslieMatrix on model-deterministic.R
# lambda.numericDistri on models-compoundDistributions.R

# Proportions for trends. Decrease includes extinct.
#' @export
trendsProp<- function(...){
  UseMethod("trendsProp")
}

#' @rdname discretePopSim
#' @export
trendsProp.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  dN<- t(diff(t(pop))) # dN =  N_t+1 - N_t
  N0<- pop[,-ncol(pop)]
  popF<- pop[,ncol(pop)] # final population
  popF[is.na(popF)]<- 0
  replicates<- nrow(pop)
  nTransitions<- sum(!is.na(dN))
  
  increase<- length(which(dN > 0)) / nTransitions
  decrease<- length(which(dN < 0)) / nTransitions
  stable<- length(which(dN == 0)) / nTransitions
  extinct<- length(which(popF == 0)) / replicates
  
  return (data.frame(increase, decrease, stable, extinct))
}

#' @rdname numericDistri
#' @export
trendsProp.numericDistri<- function(distri, N0){
  increase<- sum(distri$p[which(distri$x == (N0 + 1)):nrow(distri)])
  decrease<- sum(distri$p[1:which(distri$x == (N0 - 1))])
  stable<- distri$p[distri$x == N0]
  extinct<- distri$p[distri$x == 0]
  return(data.frame(increase, decrease, stable, extinct))
}

## Graphics ----

#' Plot population size time series of a discretePopSim simulation with replicates.
#' 
#' @rdname discretePopSim
#' @param x a discretePopSim object.
#' @param ... parameters to \code{\link[graphics]{matplot}}.
#'
#' @return
#' @export
plot.discretePopSim<- function(x, type="l", xlab="t", ylab="N", ...){
  x<- t(x)
  graphics:: matplot(x, type=type, xlab=xlab, ylab=ylab, ...)
}

#' Plot a histogram with the final population size of each replicate.
#'
#' @rdname discretePopSim
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
hist.discretePopSim<- function(x, xlab="N", ...){
  main<- as.expression(bquote("N_t=" * .(ncol(x) - 1) * " for " * .(nrow(x)) * " replicates", where=environment()))
  x<- x[,ncol(x)]
  graphics::hist(x, main=main, xlab=xlab, ...)
}


## Ideas: Lyapunov exponents
