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
#' @param raw if \code{TRUE} return the complete time serie, otherwise only t0  and tf. Useful when memory is limited.
#'
#' @return
#' @export
#'
#' @examples
discreteABMSim<- function(N0=c(N1s=5, N1b=5, N1bF=5, N2s=5, N2b=5, N2bF=5),
                          transitionsFunc=transitionABM.LH_Beh,
                          params=list(b1=2, b2=2,   broods=2, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                          tf=10, replicates=100, maxN=10000, Ntf=FALSE){
  stateName<- names(N0)
  nStates<- length(stateName)
  
  if (!Ntf){
    popABM<- array(NA_real_, dim=c(replicates, nStates, tf+1), dimnames=list(replicate=NULL, state=stateName, t=0:tf))
    popABM[,,1]<- N0
    
    for (ti in 1:tf){
      popABM[,,ti+1]<- transitionsFunc(N=popABM[,,ti], params=params)
      popABM[,,ti+1]<- apply(popABM[,,ti+1], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      if (anyNA(popABM[,,ti+1])){
        warning("NAs produced during the simulation of the discreteABMSim model.")
        if (ti < tf) popABM[,,(ti+2):(tf+1)]<- -1
        else popABM[,,tf+1]<- -1
        break
      }
      # Stop if all replicates have a class that reach maxN
      if (all(apply(popABM[,,ti+1], MARGIN=1, function(x) any(c(x == maxN, FALSE), na.rm=TRUE)))){ # FALSE in case all is NA
        popABM[,,tf+1]<- maxN
        break
      }# Stop if all replicates get extinct
      if (all(apply(popABM[,,ti+1], MARGIN=1, function(x) all(c(x <= 0, TRUE), na.rm=TRUE))) & ti < tf){ # TRUE in case all is NA
        popABM[,,(ti+2):(tf+1)]<- 0
        break
      }
    }
  }else{ # Save t0 and tf and discard intermediate timesteps
    popABM<- array(NA_real_, dim=c(replicates, nStates, 2), dimnames=list(replicate=NULL, state=stateName, t=c(0, tf)))
    popABM[,,1]<- N0
    popABM[,,2]<- N0
    
    for (ti in 1:tf){
      popABM[,,2]<- transitionsFunc(N=popABM[,,2], params=params)
      popABM[,,2]<- apply(popABM[,,2], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))

      # Stop if all replicates have a class that reach maxN
      if (all(apply(popABM[,,2], MARGIN=1, function(x) any(c(x == maxN, FALSE), na.rm=TRUE)))){ # FALSE in case all is NA
        popABM[,,2]<- maxN
        pop[,ti + 1]<- maxN
        break
      }# Stop if all replicates get extinct
      if (all(apply(popABM[,,2], MARGIN=1, function(x) all(c(x <= 0, TRUE), na.rm=TRUE)))){ # TRUE in case all is NA
        break
      }
    }
  }
  
  
  pop<- discreteABMSim2discretePopSim(popABM) # sort replicates by final size
  
  popABM<- popABM[order(pop[,ncol(pop)]),,, drop=FALSE]
  
  if (Ntf){
    popABM<- popABM[,,c(1, dim(popABM)[3]), drop=FALSE]
  }
  
  class(popABM)<- c("discreteABMSim", "array")
  return(popABM)
}


#' discreteABMSim2discretePopSim
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
  
  x<- rowSums(x[,,dim(x)[3]])
  
  graphics::hist(x, main=main, xlab=xlab, ...)
}


