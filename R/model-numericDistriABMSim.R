#' Discrete Time Agent Based Model
#' 
#' \code{numericDistriABMSim} class represent the result of numeric probability distribution simulations with
#' replicates on rows, state on columns and time on the third dimension. The class inherits from \code{array} for 
#' subsetting by replicates, state or timesteps.
#' 
#' @name numericDistriABMSim
NULL

#' numericDistriABMSim
#'
#' @rdname numericDistriABMSim
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
#' @return a \code{numericDistriABMSim} object.
#' @export
#'
#' @examples
numericDistriABMSim<- function(N0=c(N1s=0, N1b=1, N1bF=0, N2s=0, N2b=1, N2bF=0),
                          transitionsFunc=LHR:::transitionABM.LH_Beh_DIST,
                          params=list(b1=4, b2=4,   broods=1, PbF1=.4, PbF2=.4,  a1=.3,ab1=.25,sa1=.25,j1=.1,  a2=.3,ab2=.25,sa2=.20,j2=.1, AFR=1, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5),
                          tf=3, maxN=100000, Ntf=FALSE, randomizeN0=FALSE){
  # Check nStates returned by transitionsFunc (LH_behavior add subadult classes according to AFR)
  N<- transitionsFunc(N=N0, params=params)
  stateName<- names(N)
  nStates<- length(stateName)
  
  if (sum(N0 > 0) < 2) randomizeN0<- FALSE
  
  if (randomizeN0){
    N0rand<- t(N0)
    N0rand<- apply(N0rand, 2, function(x){
      selNon0<- which(x > 0)
      selNon0rand<- sample(selNon0, size=length(selNon0))
      x[selNon0]<- x[selNon0rand]
      x
    })
  }
  
  if (!Ntf){
    distriABM<- structure(as.list(rep(NA, tf+1)), names=0:tf)
    
    if (randomizeN0){
      distriABM[[1]]<- N0rand
    }else{
      distriABM[[1]]<- t(N0)
    }
    
    for (ti in 1:tf){
      distriABM[[ti+1]]<- transitionsFunc(N=distriABM[[ti]], params=params)
      
      rangeDistri<- sapply(distriABM[[ti+1]], function(x) range(x$x))
      if (max(rangeDistri) > maxN){
        warning("Population size reached maxN.")
        break
      }
      # if (min(rangeDistri) > maxN){
      #   warning("Population size reached maxN.")
      #   break
      # }
    }
  }else{ # Save t0 and tf and discard intermediate timesteps
    distriABM<- structure(list(NA, NA), names=c(0, tf))
    
    if (randomizeN0){
      distriABM[[1]]<- N0rand
      distriABM[[2]]<- N0rand
    }else{
      distriABM[[1]]<- t(N0)
      distriABM[[2]]<- t(N0)
    }

    
    for (ti in 1:tf){
      distriABM[[2]]<- transitionsFunc(N=distriABM[[2]], params=params)
      # distriABM[,, 2]<- apply(distriABM[,, 2, drop=FALSE], MARGIN=2, function(x) ifelse(x > maxN, maxN, x))
      
      rangeDistri<- sapply(distriABM[[2]], function(x) range(x$x))
      if (max(rangeDistri) > maxN){
        warning("Population size reached maxN.")
        break
      }
      # if (min(rangeDistri) > maxN){
      #   warning("Population size reached maxN.")
      #   break
      # }
    }
  }
  
  class(distriABM)<- c("numericDistriABMSim", "list")
  return(distriABM)
}


#' numericDistriABMSim2numericDistriSim
#'
#' @param distriABM 
#' @param maxN
#' @param omitClass character string containing a regular expression to exclude classes.
#'
#' @return
#' @export
#'
#' @examples
numericDistriABMSim2numericDistriSim<- function(distriABM, maxN, omitClass){
  if (!missing(omitClass)){
    distriABM<- lapply(distriABM, function(x){
          x[!grepl(omitClass, names(x))]
        })
  }
  
  distri<- lapply(distriABM, function(x){
    if (inherits(x[[1]], "numericDistri")){
      cmd<- paste(paste0("x[['", names(x), "']]"), collapse="+")
      res<- eval(parse(text=cmd))
    }else{
      res<- sum(x)
    }
    
    return(res)
  })
  
  # distri<- cleannumericDistriSim(distri)
  
  class(distri)<- c("numericDistriPopSim")

  return(distri)
}


## TODO: Graphics ----
# Plot population size time series of a numericDistriABMSim simulation with replicates.
# 
# @rdname numericDistriABMSim
# @param x a numericDistriABMSim object.
# @param ... parameters passed to \code{\link[graphics]{matplot}}.
#
# @return
# @export
#
# @examples
plot.numericDistriABMSim<- function(x, type="l", xlab="t", ylab="N", ...){
  x<- t(apply(x, MARGIN=3, rowSums, dims=1))
  graphics::matplot(x, type=type, xlab=xlab, ylab=ylab, ...)
}

# Plot a histogram with the final population size of each replicate.
# 
# @rdname numericDistriABMSim
# @param x 
# @param ... parameters passed to \code{\link[graphics]{hist}}.
#
# @return
# @export
#
# @examples
hist.numericDistriABMSim<- function(x, main, xlab="N", ...){
  if (missing(main))
    main<- as.expression(bquote("N"[t] == .(dim(x)[3] - 1) * " for " * .(dim(x)[1]) * " replicates", where=environment()))
  
  x<- x[,,dim(x)[3]]
  x<- rowSums(x, na.rm=TRUE)
  
  graphics::hist(x, main=main, xlab=xlab, ...)
}


