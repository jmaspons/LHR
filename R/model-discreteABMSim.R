#' Discrete Time Agent Based Model
#' 
#' \code{discreteABMSim} class represent the result of discrete time simulations with
#' replicates on rows, state on columns and time on the third dimension. The class inherits from \code{array} for 
#' subsetting by replicates, state or timesteps.
#' 
#' @name discreteABMSim
NULL

discreteABMSim<- function(N0=c(N1s=5, N1b=5, N1bF=5, N1j=5, N2s=5, N2b=5, N2bF=5, N2j=5),
                          transitionsFunc=transitionABM.LH_Beh,
                          params=list(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                          tf=10, replicates=10, maxN=100000){
  stateName<- names(N0)
  nStates<- length(stateName)
  
  popABM<- array(NA_real_, dim=c(replicates, nStates, tf+1), dimnames=list(replicate=NULL, state=stateName, t=0:tf))
  popABM[,1:nStates,1]<- N0
  
  for (ti in 1:tf){
    N<- popABM[,1:nStates,ti]
    popABM[,1:nStates,ti+1]<- transitionsFunc(N=N, params=params)
    # popABM[,which(popABM[,1:nStates,ti+1] > maxN),t+1]<- maxN # TODO
  }
  
  pop<- discreteABMSim2discretePopSim(popABM)
  
  popABM<- popABM[order(pop[,ncol(pop)]),,]
  class(popABM)<- c("discreteABMSim", "array")
  
  return(popABM)
}


discreteABMSim2discretePopSim<- function(popABM){
  pop<- apply(popABM, MARGIN="t", rowSums)
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")

  return(pop)
}




## Explore ABM ----
## Call
# res<- exploreABM(x0L=x0L, params=params, transitionsFunc=transitionsFunc, 
#                  tf=tf, replicates=replicates, discretePop=discretePop, finalPop=finalPop, cl=cl, ...)

exploreABM<- function(x0L=c(N1s=5, N1b=5, N1bF=5, N1j=5, N2s=5, N2b=5, N2bF=5, N2j=5),
                      params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=c(.1,.8)  ,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                      transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=100,
                      raw=TRUE, discretePop=TRUE, finalPop=TRUE, burnin=-1, maxN=100000,
                      cl=parallel::detectCores(), verbose=FALSE, ...){
  if (is.numeric(x0L)){
    x0L<- list(x0L)
  }
  
  # resStats<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=12,
  #                                 dimnames=list(scenario=paste0(rep(rownames(params), each=length(x0L)), "_N", rep(sapply(x0L, sum), times=nrow(params))),
  #                                               stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))))
  resStats<- data.frame()
  
  if (raw) resABM<- list()
  if (discretePop) resPop<- list()
  if (finalPop){
    # Ntf<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=replicates, dimnames=list(scenario=dimnames(resStats)[[1]], Nf=NULL)))
    Ntf<- data.frame()
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
    # transitions<- transitionMat(params[i,]) -> transitionsFunc
    
    if (discretePop) resPop[[i]]<- list()
    
# for (j in seq_along(x0L)){
#   N0<- x0L[[j]]
    pars<- params[i,]
    
    parallel::clusterExport(cl=cl, c("discreteABMSim", "discreteABMSim2discretePopSim", "transitionsFunc", "pars", "replicates", "tf", "maxN", "raw"), envir=environment())
    # parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
    parallel::clusterEvalQ(cl, library(LHR))

    simL<- parallel::parLapply(cl=cl, x0L, function(x){
      sim<- discreteABMSim(N0=x, transitionsFunc=transitionsFunc, params=pars, tf=tf, replicates=replicates, maxN=maxN)

      if (raw) return(sim)
      else return (discreteABMSim2discretePopSim(sim))
    })
    
    # simL<- lapply(x0L, function(x){
    #   sim<- discreteABMSim(N0=x, transitionsFunc=transitionsFunc, params=pars, tf=tf, replicates=replicates, maxN=maxN)
    #   if (raw) return(sim)
    #   else return (discreteABMSim2discretePopSim(sim))
    # })
    
    if (raw){
      popL<- lapply(simL, discreteABMSim2discretePopSim)
      resABM[[i]]<- simL
    }else{
      popL<- simL
    }
    
    # if (length(unique(sapply(simL, ncol))) > 1){warning("Discrete simulation results differs in time steps")} ## important when !is.null(dtDiscretize)
    # pop<- do.call("rbind", simL) ## TODO: Check different times for extinctions for models-discreteTime.R
    
    # Save discrete population stats
    # pop: dt= dtDiscretize
    # popDtF: dt= tf - 0
    popDtFL<- lapply(popL, function(pop){
      pop<- pop[,c(1,ncol(pop))]
      pop[is.na(pop[,2]),2]<- 0
      class(pop)<- "discretePopSim"
      return(pop)
    })
    
    stats<- lapply(popDtFL, summary)
    stats<- do.call(rbind, stats)
    
    N0<- sapply(x0L, sum, na.rm=TRUE)
    
    tmp<- data.frame(scenario=rownames(params)[i], N0=N0, stats, stringsAsFactors=FALSE) ## TODO check params column. It's always 1!!
    if (verbose) print(tmp, row.names=FALSE)
    #       if (!is.null(dtDiscretize)) tmp<- rbind(tmp, summary(pop))
    #       names(tmp)<- gsub(".1", ".dt", names(tmp))
    # if(paste0(rownames(params)[i], "_N", N0) != rownames(resStats)[k]) stop("Incorrect rownames")
    resStats<- rbind(resStats, tmp)
    

    
    if (discretePop) resPop[[i]]<- popL
    if (finalPop){
      popTfL<- lapply(popL, function(pop){
        pop[,ncol(pop)]
        pop[is.na(pop)]<- 0
        sort(pop)
      })
      Ntf<- rbind(Ntf, do.call(rbind, popTfL))
    }
    k<- k +1
# }# End N0 loop
    
    if (discretePop) names(resPop[[i]])<- paste0("N", N0)
  }# End parameters loop
  
  #   names(resStats)<- rownames(params)
  
  res<- list(stats=resStats)
  res<- c(res, list(params=params))
  res<- c(res, list(simParams=list(x0L=x0L, tf=tf, replicates=replicates, burnin=burnin)))
  
  if (raw){
    res<- c(res, list(discreteABMSim=resABM))
  }
  if (discretePop){
    names(resPop)<- rownames(params)
    res<- c(res, list(discretePopSim=resPop))
  }
  if (finalPop){
    res<- c(res, list(Ntf=Ntf))
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  class(res)<- "exploreABMSim"
  return (res)
}


