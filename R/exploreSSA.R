#' Class ssa
#' 
#' @name ssa
#' @exportClass ssa
NULL

## Explore the deterministic part of the model ----
# For inspiration chapter 3: Tuljapurkar, Shripad, and Hal Caswell, eds. Structured-population models in marine, terrestrial, and freshwater systems. Vol. 18. Springer Science & Business Media, 2012.
# Simulate the model with discrete time and continuous population size
# Useful to find the stable state structure and deterministic growth rate
# normalize=TRUE work with the proportions of states instead of the population size.
#   It converges faster and allows to calculated the maximum grow rates for models with densodependence,
#   assuming that the carring capacity >> 1
ssa.deterministic<- function(init.values=rep(1, nrow(transitions)), transitions, rateFunc, params, maxTf=500, normalize=TRUE){
  res<- matrix(NA_real_, nrow=maxTf, ncol=9, dimnames=list(NULL, state=c("time", "N1s",  "N1b",  "N1bF", "N1j",  "N2s",  "N2b",  "N2bF", "N2j")))

  res[1,]<- c(1, init.values / sum(init.values))
  converge<- FALSE
  i<- 1
  while (i < maxTf & !converge){
    i<- i + 1
    freqs<- rateFunc(res[i-1,], params=params)
    if (normalize) freqs<- freqs / sum(freqs)
    delta<- t(freqs * t(transitions))
    delta<- rowSums(delta)
    
    N<- res[i-1,-1] + delta
    N[which(N < 0)]<- 0
    if (sum(N) == 0){
      res[i,]<- c(i, N)
      break
    }
    
    if (normalize) N<- N / sum(N)
    res[i,]<- c(i, N)
    if (sum(diff(res[c(i-1, i), -1])) == 0) converge<- TRUE
  }
  
  stableStructure<- res[i,-1]
  lambda<- sum(delta) + 1 # Nt+1 / Nt
  r<- sum(delta) # (Nt+1 - Nt) / Nt
  
  return (list(stableStructure=stableStructure, lambda=lambda, r=r))
}


## Explore SSA ----
exploreSSA<- function(x0L, params, transitionMat, rateFunc, maxTf=10, replicates=100,
                      discretePop=FALSE, finalPop=TRUE, burnin=-1, dtDiscretize=NULL, 
                      cl=parallel::detectCores(),
                      verbose=FALSE, ...){
  if (is.numeric(x0L)){
    x0L<- list(x0L)
  }

  resStats<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=12,
                                  dimnames=list(scenario=paste0(rep(rownames(params), each=length(x0L)), "_N", rep(sapply(x0L, sum), times=nrow(params))),
                                                stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))))
  if (discretePop) resPop<- list()
  if (finalPop){
    Ntf<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=replicates, dimnames=list(scenario=dimnames(resStats)[[1]], Nf=NULL)))
  }
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }
  
  k<- 1
  for (i in 1:nrow(params)){
    message(i, "/", nrow(params), "\t", rownames(params)[i])
    transitions<- transitionMat(params[i,])
    
    if (discretePop) resPop[[i]]<- list()

    for (j in seq_along(x0L)){
      N0<- x0L[[j]]
      pars<- params[i,]

      parallel::clusterExport(cl=cl, c("N0", "transitions", "rateFunc", "pars", "maxTf", "dtDiscretize", "burnin"), envir=environment())
      parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
      parallel::clusterEvalQ(cl, library(LHR))
      
      simL<- parallel::parLapply(cl=cl, 1:replicates, function(x){
        sim<- adaptivetau::ssa.adaptivetau(init.values=N0, transitions=transitions, rateFunc=rateFunc, params=pars, tf=maxTf, ...)
        sim<- discretizePopSim(sim, dt=dtDiscretize, burnin=burnin)
        return (sim)
      })

      if (length(unique(sapply(simL, ncol))) > 1){warning("Discrete simulation results differs in time steps")} ## important when !is.null(dtDiscretize)
      pop<- do.call("rbind", simL) ## TODO: Check different times for extinctions for models-discreteTime.R

      # Save discrete population stats
      # pop: dt= dtDiscretize
      # popDtF: dt= maxTf - 0
      popDtF<- pop[,c(1,ncol(pop))]
      popDtF[is.na(popDtF[,2]),2]<- 0
      tmp<- data.frame(scenario=rownames(params)[i], N0=sum(x0L[[j]]), summary(popDtF), stringsAsFactors=FALSE) ## TODO check params column. It's always 1!!
#       if (!is.null(dtDiscretize)) tmp<- rbind(tmp, summary(pop))
#       names(tmp)<- gsub(".1", ".dt", names(tmp))
if(paste0(rownames(params)[i], "_N", sum(x0L[[j]])) != rownames(resStats)[k]) stop("Incorrect rownames")
      resStats[k,]<- tmp

      if (verbose) print(tmp, row.names=FALSE)

      if (discretePop) resPop[[i]][[j]]<- pop
      if (finalPop){
        popTf<- pop[,ncol(pop)]
        popTf[is.na(popTf)]<- 0

        Ntf[k,]<- sort(popTf)
      }
      k<- k +1
    }# End N0 loop

    if (discretePop) names(resPop[[i]])<- paste0("N", sapply(x0L, sum))
  }# End parameters loop

#   names(resStats)<- rownames(params)

  res<- list(stats=resStats)
  res<- c(res, list(params=params))
  res<- c(res, list(simParams=list(x0L=x0L, tf=maxTf, replicates=replicates, burnin=burnin, dtDiscretize=dtDiscretize)))

  if (discretePop){
    names(resPop)<- rownames(params)
    res<- c(res, list(pop=resPop))
  }
  if (finalPop){
    res<- c(res, list(Ntf=Ntf))
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  class(res)<- "ssa"
  return (res)
}


exploreSSA.deterministic<- function(x0=rep(1, nrow(transitionMat())), params, transitionMat, rateFunc, maxTf=500, normalize=TRUE,
                                    cl=parallel::detectCores(), ...){
  params<- split(params, rownames(params))
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }

  parallel::clusterExport(cl=cl, c("x0", "transitionMat", "rateFunc", "maxTf", "normalize"), envir=environment())
  parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
  parallel::clusterEvalQ(cl, library(LHR))
  
  res<- parallel::parLapply(cl=cl, params, function(x){
    transitions<- transitionMat(x)
    sim<- ssa.deterministic(init.values=x0, transitions=transitions, rateFunc=rateFunc, params=x, maxTf=maxTf, normalize=normalize)
    sim<- c(lambda=sim$lambda, r=sim$r, sim$stableStructure)
    return (sim)
  })
  
  res<- as.data.frame(do.call("rbind", res))
  
  if (numericCL) parallel::stopCluster(cl)
  
  return (res)
}


#' Discretize the result of a ssa simulations from adaptivetau
#' @param sim matrix(ncol=1:tf, nrow=nStates, dimnames=list(time=NULL, state=c("time", states))) from ssa (e.g. models-IBM-ssa_LH_behavior.R)
#' @param dt 
#' @param burnin 
#' @param keepStates 
#'
#' @return a \code{\link{discretePopSim}} object.
#' @export
#'
#' @examples
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
    class(pop)<- c("discretePopSim")
    return(pop)
  }
}
