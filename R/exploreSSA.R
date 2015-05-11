setOldClass("ssa")

## Explore SSA
exploreSSA<- function(x0L, params, transitionMat, rateFunc, tf=10, replicates=100,
                      discretePop=FALSE, finalPop=TRUE, burnin=-1, dtDiscretize=NULL, cores=1, mc.preschedule=TRUE, ...){
  if (is.numeric(x0L)){
    x0L<- list(x0L)
  }

  resStats<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=12,
                                          dimnames=list(scenario=paste0(rownames(params), "_N", rep(sapply(x0L, sum), times=nrow(params))),
                                                                        stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))))
#                                                         , "increase.dtF", "decrease.dtF", "stable.dtF", "extinct.dtF", "GR.dtF", "meanR.dtF", "varR.dtF", "GL.dtF", "meanL.dtF", "varL.dtF"))))
  if (discretePop) resPop<- list()
  if (finalPop) fPop<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=replicates, dimnames=list(scenari=paste0(rownames(params), "_N", rep(sapply(x0L, sum), times=nrow(params)), Nf=NULL))))
  k<- 1
  for (i in 1:nrow(params)){
    cat(i, "/", nrow(params), rownames(params)[i], "\n")
    transitions<- transitionMat(params[i,])
    
    if (discretePop) resPop[[i]]<- list()
#     if (finalPop) fPop[[i]]<- as.data.frame(matrix(nrow=length(x0L), ncol=replicates, dimnames=list(N0=sapply(x0L, sum), Nf=NULL)))

    for (j in seq_along(x0L)){
      RNGkind("L'Ecuyer-CMRG") # ?mcparallel > Random numbers
      mc.reset.stream()
      simL<- mclapply(1:replicates, function(x){
        sim<- ssa.adaptivetau(init.values=x0L[[j]], transitions=transitions, rateFunc=rateFunc, params=params[i,], tf=tf)# TODO, ...)
        sim<- discretizePopSim(sim, dt=dtDiscretize, burnin=burnin)
        return (sim)
      }, mc.cores=cores, mc.set.seed=TRUE, mc.preschedule=mc.preschedule)

      if (length(unique(sapply(simL, ncol))) > 1){warning("Discrete simulation results differs in time steps")} ## important when !is.null(dtDiscretize)
      pop<- do.call("rbind", simL) ## TODO: Check different times for extinctions for models-discreteTime.R

      # Save discrete population stats
      # pop: dt= dtDiscretize
      # popDtF: dt= tf - 0
      popDtF<- pop[,c(1,ncol(pop))]
      popDtF[is.na(popDtF[,2]),2]<- 0
      tmp<- data.frame(scenario=rownames(params)[i], N0=sum(x0L[[j]]), summary(popDtF), stringsAsFactors=FALSE) ## TODO check params column. It's always 1!!
#       if (!is.null(dtDiscretize)) tmp<- rbind(tmp, summary(pop))
#       names(tmp)<- gsub(".1", ".dt", names(tmp))
      resStats[k,]<- tmp

      print(tmp, row.names=FALSE)

      if (discretePop) resPop[[i]][[j]]<- pop
      if (finalPop) fPop[k,]<- sort(pop[,ncol(pop)])
      k<- k +1
    }# End N0 loop

    if (discretePop) names(resPop[[i]])<- paste0("N", sapply(x0L, sum))
  }# End parameters loop

#   names(resStats)<- rownames(params)

  res<- list(stats=resStats)
  if (discretePop){
    names(resPop)<- rownames(params)
    res<- c(res, list(pop=resPop))
  }
  if (finalPop){
    res<- c(res, list(popTf=fPop))
  }
  res<- c(res, list(params=params))
  res<- c(res, list(simParams=list(x0L=x0L, tf=tf, replicates=replicates, burnin=burnin, dtDiscretize=dtDiscretize)))
  class(res)<- "ssa"
  return (res)
}


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
