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
#' @param Ntf if \code{FALSE} return the complete time serie, otherwise only t0  and tf. Useful when memory is limited.
#' @param randomizeN0 randomly exchange N0 for classes > 0. Useful for findN0 whith non balanced N0.
#'
#' @return a \code{discreteABMSim} object.
#' @export
#'
#' @examples
discreteABMSim<- function(N0=c(N1s=5, N1b=5, N1bF=5, N2s=5, N2b=5, N2bF=5),
                          transitionsFunc=transitionABM.LH_Beh,
                          params=list(b1=2, b2=2,   broods=2, PbF1=.4, PbF2=.4,  a1=.3,ab1=.25,sa1=.25,j1=.1,  a2=.3,ab2=.25,sa2=.20,j2=.1, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5),
                          tf=10, replicates=100, maxN=100000, Ntf=FALSE, randomizeN0=FALSE){
  # Check nStates returned by transitionsFunc (LH_behavior add subadult classes according to AFR)
  N<- transitionsFunc(N=rbind(N0, N0), params=params)
  stateName<- colnames(N)
  nStates<- length(stateName)
  
  if (length(which(N0 > 0)) < 2) randomizeN0<- FALSE
  
  if (!Ntf){
    popABM<- array(0, dim=c(replicates, nStates, tf+1), dimnames=list(replicate=NULL, state=stateName, t=0:tf))
    
    if (randomizeN0){
      N0rand<- replicate(replicates, N0)
      N0rand<- apply(N0rand, 2, function(x){
        selNon0<- which(x > 0)
        selNon0rand<- sample(selNon0, size=length(selNon0))
        x[selNon0]<- x[selNon0rand]
        x
      })
      popABM[, names(N0), 1]<- N0rand
    }else{
      popABM[, names(N0), 1]<- t(replicate(replicates, N0))
    }
    
    for (ti in 1:tf){
      popABM[,, ti+1]<- transitionsFunc(N=popABM[,, ti], params=params)
      popABM[,, ti+1]<- apply(popABM[,, ti+1], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      if (anyNA(popABM[,, ti+1])){
        warning("NAs produced during the simulation of the discreteABMSim model.")
        if (ti < tf) popABM[,, (ti+2):(tf+1)]<- -1
        else popABM[,, tf+1]<- -1
        break
      }
      
      if (ti %% 10 == 0 & ti + 1 < tf){ # check stop conditions every 10 time steps
        # Stop if all replicates have a class that reach maxN. TODO: check if the optimization is worth it benchmark.Rmd
        if (all(apply(popABM[,, ti+1], MARGIN=1, function(x) any(c(x == maxN, FALSE), na.rm=TRUE)))){ # FALSE in case all is NA
          popABM[,, ti+1]<- maxN
          popABM[,, (ti+2):(tf+1)]<- NA # remove 0. If maxN is not stable, the transitions are cosidered valid at maxNNA(pop)
          break
        }# Stop if all replicates get extinct TODO: check if the optimization is worth in benchmark.Rmd
        if (all(apply(popABM[,, ti+1], MARGIN=1, function(x) all(c(x <= 0, TRUE), na.rm=TRUE))) & ti < tf){ # TRUE in case all is NA
          break
        }
      }
    }
  }else{ # Save t0 and tf and discard intermediate timesteps
    popABM<- array(0, dim=c(replicates, nStates, 2), dimnames=list(replicate=NULL, state=stateName, t=c(0, tf)))
    
    if (randomizeN0){
      N0rand<- replicate(replicates, N0)
      N0rand<- apply(N0rand, 2, function(x){
        selNon0<- which(x > 0)
        selNon0rand<- sample(selNon0, size=length(selNon0))
        x[selNon0]<- x[selNon0rand]
        x
      })
      popABM[, names(N0), 1]<- N0rand
      popABM[, names(N0), 2]<- N0rand
    }else{
      popABM[, names(N0), 1]<- t(replicate(replicates, N0))
      popABM[, names(N0), 2]<- t(replicate(replicates, N0))
    }

    
    for (ti in 1:tf){
      popABM[,, 2]<- transitionsFunc(N=popABM[,, 2], params=params)
      popABM[,, 2]<- apply(popABM[,, 2], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      if (ti %% 10 == 0 & ti < tf){ # check stop conditions every 10 time steps
        # Stop if all replicates have a class that reach maxN
        if (all(apply(popABM[,, 2], MARGIN=1, function(x) any(c(x == maxN, FALSE), na.rm=TRUE)))){ # FALSE in case all is NA
          break
        }# Stop if all replicates get extinct
        if (all(apply(popABM[,, 2], MARGIN=1, function(x) all(c(x <= 0, TRUE), na.rm=TRUE)))){ # TRUE in case all is NA
          popABM[,, 2]<- 0
          break
        }
      }
    }
  }
  
  
  pop<- discreteABMSim2discretePopSim(popABM) # sort replicates by final size
  
  popABM<- popABM[order(pop[, ncol(pop)], na.last=TRUE),,, drop=FALSE]
  
  class(popABM)<- c("discreteABMSim", "array")
  return(popABM)
}


#' discreteABMSim2discretePopSim
#'
#' @param popABM 
#' @param maxN
#' @param omitClass character string containing a regular expression to exclude classes.
#'
#' @return
#' @export
#'
#' @examples
discreteABMSim2discretePopSim<- function(popABM, maxN, omitClass){
  if (missing(omitClass)){
    pop<- apply(popABM, MARGIN=3, rowSums)
  }else{
    pop<- apply(popABM[, !grepl(omitClass, colnames(popABM)),, drop=FALSE], MARGIN=3, rowSums)
  }
  
  
  pop<- cleanDiscretePopSim(pop)
  
  class(pop)<- c("discretePopSim", "matrix")

  return(pop)
}


## Graphics ----
#' Plot population size time series of a discreteABMSim simulation with replicates.
#' 
#' @rdname discreteABMSim
#' @param x a discreteABMSim object.
#' @param ... parameters passed to \code{\link[graphics]{matplot}}.
#'
#' @return
#' @export
#'
#' @examples
plot.discreteABMSim<- function(x, type="l", xlab="t", ylab="N", ...){
  x<- t(apply(x, MARGIN=3, rowSums, dims=1))
  graphics::matplot(x, type=type, xlab=xlab, ylab=ylab, ...)
}

#' Plot a histogram with the final population size of each replicate.
#' 
#' @rdname discreteABMSim
#' @param x 
#' @param ... parameters passed to \code{\link[graphics]{hist}}.
#'
#' @return
#' @export
#'
#' @examples
hist.discreteABMSim<- function(x, main, xlab="N", ...){
  if (missing(main))
    main<- as.expression(bquote("N"[t] == .(dim(x)[3] - 1) * " for " * .(dim(x)[1]) * " replicates", where=environment()))
  
  x<- x[,,dim(x)[3]]
  x<- rowSums(x, na.rm=TRUE)
  
  graphics::hist(x, main=main, xlab=xlab, ...)
}


