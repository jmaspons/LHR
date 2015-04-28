# sim: matrix(ncol=1:tf, nrow=nStates, dimnames=list(time=NULL, state=c("time", states))) from ssa (e.g. models-IBM-ssa_LH_behavior.R)
discretizePopSim<- function(sim, dt=NULL, burnin=-1, keepStates=FALSE){
  extinctionTime<- sim[match(TRUE, rowSums(sim[,-1]) == 0), "time"]
  sim<- sim[sim[,"time"] > burnin,]
  
  if (nrow(sim) == 0) return (NA) # if population goes extinct before burnin return NA

  if (is.null(dt)){ # keep first an last time only
    sim<- sim[c(1, nrow(sim)),]
    period<- ordered(c(sim[1,"time"], sim[nrow(sim), "time"]))
    steps<- c(sim[1,"time"], sim[nrow(sim), "time"])
  }else{
    steps<- seq(sim[1,"time"], sim[nrow(sim), "time"], by=dt)
    if (length(steps) < 2){
      sim<- sim[c(1, nrow(sim)),] # Pick the first and last values
      period<- ordered(c(sim[1,"time"], sim[nrow(sim), "time"]))
      steps<- c(sim[1,"time"], sim[nrow(sim), "time"])
    }else{
      period<- cut(sim[,"time"], steps, include.lowest=TRUE, right=FALSE, ordered_result=FALSE)
    }
  }

  selSteps<- sapply(levels(period), match, period)
  sim<- sim[selSteps,]
  sim[,"time"]<- steps[1:nrow(sim)]

  if (!is.na(extinctionTime)){
    extinctionStep<- match(TRUE, extinctionTime < sim[,"time"]) - 1 ## WARNING < or <= ??
    if (!is.na(extinctionStep)){ # No extinction on the last timestep
      if (sum(sim[extinctionStep, -1]) > 0 & extinctionStep < nrow(sim)){
      # If population goes extinct and the event doesn't appear on the selected timesteps (e.g. the extinction occurs after the first event of the interval)
      #  add population = 0 after the last interval with population > 0
        extinctionStep<- extinctionStep + 1
        sim[extinctionStep, -1]<- 0
      }
      if (extinctionStep < nrow(sim)){
        sim[(extinctionStep+1):nrow(sim), -1]<- NA # fill with NA after the extinction
      }
    }
    noEvents<- which(is.na(sim[,2]) & sim[,"time"] < extinctionTime) # time steps with no events
  }else{
    noEvents<- which(is.na(sim[,2]) & sim[,"time"]) # time steps with no events
  }

  # fill intervals with no events with the same value as the former step
  for (i in noEvents){
    sim[i,-1]<- sim[i-1,-1]
  }
#   noEventsOK<- noEvents[c(1, which(diff(noEvents) > 1) + 1)] # time steps with a valid time step just before
#   noEvents2<- noEvents[which(diff(noEvents) == 1) + 1] # No events just before
#   diff(noEvents2)
#   sim[noEventsOK,-1]<- sim[noEventsOK - 1,-1]

  if (keepStates){
    return (sim)
  }else{
    pop<- matrix(rowSums(sim[,-1]), nrow=1, dimnames=list(replicate=NULL, t=sim[,1]))
    pop<- as.data.frame(pop)
    class(pop)<- c("discretePopSim", "data.frame")
    return(pop)
  }
}

# pop: matrix(nrow=replicates, ncol=tf+1, dimnames=list(replicate=NULL, t=0:tf)) from models-discreteTime.R
# pop<- mSurvBV.t(broods=2, clutch=2, nestFail=.5, juvSurv=.5, adultSurv=.7, N0=20, replicates=100, tf=10)
summary.discretePopSim<- function(pop, dt=1){
  R<- unlist(r(pop, dt=dt)) # intrinsic growth rate
  L<- unlist(lambda(pop, dt=dt)) # lambda
  meanR<- mean(R, na.rm=TRUE)
  varR<- var(R, na.rm=TRUE)
  meanL<- mean(L, na.rm=TRUE)
  varL<- var(L, na.rm=TRUE)
  GR<- G(meanR, varR)
  GL<- G(meanL, varL)
  trends<- trendsProp(pop)
  res<- data.frame(trends, GR, meanR, varR, GL, meanL, varL)
  return(res)
}

r<- function(...){
  UseMethod("r")  
}

r.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  dN<- t(diff(t(pop))) # dN =  N_t+1 - N_t
  N0<- pop[-ncol(pop)]
  r<- (dN / dt) / N0
#   names(r)<- colnames(dN) # otherwise it takes the names from N0 (0:(tf-1) instead of 1:tf as does lambda function
  return (r) # intrinsic grow rate (r = dN / dt / N)
}

lambda<- function(...){
  UseMethod("lambda")  
}

lambda.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  lambda<- pop[,-1] / pop[,-ncol(pop)] # lambda = Nt+1 / Nt
  return (lambda)
}

# lambda.leslieMatrix on model-deterministic.R

# Proportions for trends. Decrease includes extinct.
trendsProp<- function(...){
  UseMethod("trendsProp")
}

trendsProp.discretePopSim<- function(pop, dt=1){
  sampleT<- seq(1, ncol(pop), by=dt)
  if (length(sampleT) < 2) {warning("length(sampleT) < 2")}
  pop<- pop[, sampleT]
  dN<- t(diff(t(pop))) # dN =  N_t+1 - N_t
  N0<- pop[-ncol(pop)]
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

trendsProp.numericDistri<- function(distri, N0){
  increase<- sum(distri$p[which(distri$x == (N0 + 1)):nrow(distri)])
  decrease<- sum(distri$p[1:which(distri$x == (N0 - 1))])
  stable<- distri$p[distri$x == N0]
  extinct<- distri$p[distri$x == 0]
  return(data.frame(increase, decrease, stable, extinct))
}

## Ideas: Lyapunov exponents
