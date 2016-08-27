#' Discrete Time Agent Based Model
#' 
#' \code{discreteABMSim} class represent the result of discrete time simulations with
#' replicates on rows, state on columns and time on the third dimension. The class inherits from \code{array} for 
#' subsetting by replicates, state or timesteps.
#' 
#' @name discreteABMSim
NULL

#' discreteABMSim
#'
#' @rdname discreteABMSim
#' 
#' @param N0 
#' @param transitionsFunc 
#' @param params 
#' @param tf 
#' @param replicates 
#' @param maxN 
#'
#' @return
#' @export
#'
#' @examples
discreteABMSim<- function(N0=c(N1s=5, N1b=5, N1bF=5, N2s=5, N2b=5, N2bF=5),
                          transitionsFunc=transitionABM.LH_Beh,
                          params=list(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                          tf=10, replicates=100, maxN=10000, raw=TRUE){
  stateName<- names(N0)
  nStates<- length(stateName)
  
  if (raw){
    popABM<- array(NA_real_, dim=c(replicates, nStates, tf+1), dimnames=list(replicate=NULL, state=stateName, t=0:tf))
    popABM[,,1]<- N0
    
    for (ti in 1:tf){
      popABM[,,ti+1]<- transitionsFunc(N=popABM[,,ti], params=params)
      popABM[,,ti+1]<- apply(popABM[,,ti+1], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      # Stop if all replicates have a class that reach maxN
      if (all(apply(popABM[,,ti+1], MARGIN=1, function(x) any(x == maxN)))){
        popABM[,,tf+1]<- maxN
        break
      }# Stop if all replicates get extinct
      if (all(apply(popABM[,,ti+1], MARGIN=1, function(x) all(x <= 0))) & ti < tf){
        popABM[,,(ti+2):(tf+1)]<- 0
        break
      }
    }
  }else{
    popABM<- array(NA_real_, dim=c(replicates, nStates, 2), dimnames=list(replicate=NULL, state=stateName, t=c(0, tf)))
    popABM[,,1]<- N0
    popABM[,,2]<- N0
    
    for (ti in 1:tf){
      popABM[,,2]<- transitionsFunc(N=popABM[,,2], params=params)
      popABM[,,2]<- apply(popABM[,,2], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      # Stop if all replicates have a class that reach maxN
      if (all(apply(popABM[,,2], MARGIN=1, function(x) any(x == maxN)))){
        popABM[,,2]<- maxN
        break
      }# Stop if all replicates get extinct
      if (all(apply(popABM[,,2], MARGIN=1, function(x) all(x <= 0)))){
        break
      }
    }
  }
  
  pop<- discreteABMSim2discretePopSim(popABM)
  
  popABM<- popABM[order(pop[,ncol(pop)]),,, drop=FALSE]
  class(popABM)<- c("discreteABMSim", "array")
  
  return(popABM)
}


#' Title
#'
#' @param popABM 
#'
#' @return
#' @export
#'
#' @examples
discreteABMSim2discretePopSim<- function(popABM, omitJuv=FALSE){
  if (omitJuv){
    pop<- apply(popABM[, !grepl("j", colnames(popABM)),], MARGIN=3, rowSums)
  }else{
    pop<- apply(popABM, MARGIN=3, rowSums)
  }
  
  pop<- pop[order(pop[,ncol(pop)]),, drop=FALSE]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")

  return(pop)
}



## Call
# res<- exploreABM(x0L=x0L, params=params, transitionsFunc=transitionsFunc, 
#                  tf=tf, replicates=replicates, discretePop=discretePop, finalPop=finalPop, cl=cl, ...)

exploreABM<- function(x0L=c(N1s=5, N1b=5, N1bF=5, N2s=5, N2b=5, N2bF=5),
                      params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=c(.1,.8)  ,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                      transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=100,
                      raw=TRUE, discretePop=TRUE, finalPop=TRUE, burnin=-1, maxN=100000,
                      cl=parallel::detectCores(), verbose=FALSE, ...){
  if (is.numeric(x0L)){
    x0L<- list(x0L)
  }
  
  N0<- sapply(x0L, sum, na.rm=TRUE)
  
  resStats<- as.data.frame(matrix(nrow=length(N0) * nrow(params), ncol=12,
                                  dimnames=list(scenario=paste0(rep(rownames(params), each=length(x0L)), "_N", rep(sapply(x0L, sum), times=nrow(params))),
                                          stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))))
  # resStats<- data.frame()
  
  if (raw) resABM<- list()
  if (discretePop) resPop<- list()
  if (finalPop){
    Ntf<- as.data.frame(matrix(nrow=length(x0L) * nrow(params), ncol=replicates, dimnames=list(scenario=dimnames(resStats)[[1]], Nf=NULL)))
    # Ntf<- data.frame()
    N_ABMtf<- array(NA_real_, dim=c(replicates, length(x0L[[1]]), length(x0L) * nrow(params)),
                    dimnames=list(replicates=NULL, state=names(x0L[[1]]), scenario=dimnames(resStats)[[1]]))
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

    
    if (discretePop) resPop[[i]]<- list()
    
    pars<- params[i,]
    
    parallel::clusterExport(cl=cl, c("discreteABMSim", "discreteABMSim2discretePopSim", "transitionsFunc", "pars", "replicates", "tf", "maxN", "raw"), envir=environment())
    # parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
    parallel::clusterEvalQ(cl, library(LHR))

    simL<- parallel::parLapply(cl=cl, x0L, function(x){
      sim<- discreteABMSim(N0=x, transitionsFunc=transitionsFunc, params=pars, tf=tf, replicates=replicates, maxN=maxN, raw=raw)

      return(sim)
    })
    
    # simL<- lapply(x0L, function(x){
    #   discreteABMSim(N0=x, transitionsFunc=transitionsFunc, params=pars, tf=tf, replicates=replicates, maxN=maxN, raw=raw)
    # })
    
    names(simL)<- paste0("N", N0)
    popL<- lapply(simL, discreteABMSim2discretePopSim)
    
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
    
    
    tmp<- data.frame(scenario=rownames(params)[i], N0=N0, stats, stringsAsFactors=FALSE) ## TODO check params column. It's always 1!!
    if (verbose) print(tmp, row.names=FALSE)
    #       if (!is.null(dtDiscretize)) tmp<- rbind(tmp, summary(pop))
    #       names(tmp)<- gsub(".1", ".dt", names(tmp))
    # if(paste0(rownames(params)[i], "_N", N0) != rownames(resStats)[k]) stop("Incorrect rownames")
    resStats[k:(k + nrow(tmp) - 1),]<- tmp
    
    if (raw){
      resABM[[i]]<- simL
    }
    
    if (finalPop){
      # discreteABMSim
      popABMTfL<- lapply(simL, function(x){
        x<- x[,,dim(x)[3]]
      })

      for (n in seq_along(N0)){
        N_ABMtf[,,k+n-1]<- popABMTfL[[n]]
      }
      
      # discretePopSim
      popTfL<- lapply(popL, function(pop){
        pop<- pop[,ncol(pop)]
        pop[is.na(pop)]<- 0
        sort(pop)
      })
      
      popTfL<- do.call(rbind, popTfL)
      
      Ntf[k:(k + nrow(popTfL) - 1),]<- popTfL
    }
    
    if (discretePop){
      resPop[[i]]<- popL
    }

    k<- k + length(N0)
  }# End parameters loop
  
  
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
    res<- c(res, list(N_ABMtf=N_ABMtf))
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  class(res)<- "exploreABMSim"
  return (res)
}


