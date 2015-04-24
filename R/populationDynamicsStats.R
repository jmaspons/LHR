# pop: matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf)) from models-discreteTime.R
# pop<- mSurvBV.t(broods, clutch, nestFail, juvSurv, adultSurv, N0,replicates, tf)
summary.discretePopSim<- function(pop, dt=1){
  R<- r(pop, dt=dt) # intrinsic growth rate
  L<- lambda(pop, dt=dt) # lambda
  meanR<- mean(R)
  varR<- var(as.numeric(R))
  meanL<- mean(L)
  varL<- var(as.numeric(L))
  GR<- G(meanR, varR)
  GL<- G(meanL, varL)
  trends<- trendsProp(pop)
  res<- data.frame(trends, GR, meanR, varR, GL, meanL, varL)
  return(res)
}

r<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  dN<- t(diff(t(pop))) # dN =  N_t+1 - N_t
  return ((dN / dt) / pop[,-ncol(pop)]) # intrinsic grow rate (r = dN / dt / N)
}

lambda<- function(...){
  UseMethod("lambda")  
}

lambda.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  return (pop[,-1] / pop[,-ncol(pop)]) # lambda = Nt+1 / Nt
}

# Proportions for trends. Decrease includes extinct.
trendsProp<- function(...){
  UseMethod("trendsProp")
}

trendsProp.discretePopSim<- function(pop){
  popF<- pop[,ncol(pop)] # final population
  N0<- pop[,1]
  replicates<- nrow(pop)
  increase<- length(which(popF > N0)) / replicates
  decrease<- length(which(popF < N0)) / replicates
  stable<- length(which(popF == N0)) / replicates
  extinct<- length(which(popF == 0)) / replicates
  return (data.frame(increase, decrease, stable, extinct))
}

trendsProp.numericDistri<- function(distri, N0){
  increase<- sum(distri$p[which(distri$x == (N0 + 1)):nrow(distri)])
  decrease<- sum(distri$p[1:which(distri$x == (N0 - 1))])
  stable<- distri$p[distri$x == N0]
  extinct<- distri$p[distri$x == 0]
  return(data.frame(increase, decrease, stable, extinct))
}

## Ideas: Lyapunov exponents
