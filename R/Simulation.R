#' @include aaa-classes.R
NULL

## Constructor ----
#' @rdname Sim
#'
#' @param params 
#'
#' @return a \code{Sim} object.
#' @examples Sim()
#' @export
setGeneric("Sim", function(params, type="discretePopSim") standardGeneric("Sim"))

setMethod("Sim",
          signature(params="missing", type="ANY"),
          function(params, type="discretePopSim"){
            sim<- switch(type,
                         discretePopSim=Sim.discretePopSim(),
                         numericDistri=Sim.numericDistri(),
                         ABM=Sim.ABM(),
                         ssa=Sim.ssa())
            Sim.discretePopSim()
          }
)

setMethod("Sim",
          signature(params="Model", type="missing"),
          function(params){
            simPars<- params@sim@params

            sim<- switch(class(params@sim),
                   Sim.discretePopSim=Sim.discretePopSim(params=simPars),
                   Sim.numericDistri=Sim.numericDistri(params=simPars),
                   Sim.ABM=Sim.ABM(params=simPars),
                   Sim.ssa=Sim.ssa(params=simPars))

            return (sim)
          }
)

## Subclasses
## Sim.discretePopSim Class ----

#' @rdname Sim.discretePopSim
#'
#' @param N0 
#' @param envVar 
#' @param sexRatio 
#' @param matingSystem 
#' @param tf 
#' @param replicates 
#' @param raw 
#' @param Ntf 
#'
#' @return a \code{Sim.discretePopSim} object.
#' @examples Sim.discretePopSim()
#' 
#' @export
setGeneric("Sim.discretePopSim", function(params, N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE),
                                          sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                                          tf=10, replicates=15, raw=TRUE, Ntf=TRUE) standardGeneric("Sim.discretePopSim"))

setMethod("Sim.discretePopSim",
          signature(params="list", N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY", tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY"),
          function(params){
            sim<- new("Sim.discretePopSim", params=params)

            return (sim)
          }
)

setMethod("Sim.discretePopSim",
          signature(params="missing", N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY",
                    tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   tf=10, replicates=15, raw=TRUE, Ntf=TRUE){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }
            
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf)
            sim<- Sim.discretePopSim(params=params)
            
            return (sim)
          }
)


## Sim.numericDistri Class ----
#TODO: sexRatio

#' @rdname Sim.numericDistri
#'
#' @param N0 
#' @param envVar 
#' @param sexRatio 
#' @param matingSystem 
#' @param tf 
#' @param raw 
#'
#' @return a \code{Sim.numericDistri} object.
#' @examples Sim.numericDistri()
#' @export
setGeneric("Sim.numericDistri", function(params, N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), tf=1, raw=TRUE) standardGeneric("Sim.numericDistri"))

setMethod("Sim.numericDistri",
          signature(params="list", N0="missing", envVar="missing", sexRatio="missing", matingSystem="missing", tf="missing", raw="missing"),
          function(params){
            sim<- new("Sim.numericDistri", params=params)
            
            return (sim)
          }
)

setMethod("Sim.numericDistri",
          signature(params="missing", N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY", tf="ANY", raw="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), tf=1, raw=TRUE){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }
            
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, raw=raw)
            sim<- Sim.numericDistri(params=params)
            
            return (sim)
          }
)


## Sim.ABM Class ----

#' @rdname Sim.ABM
#'
#' @param N0 
#' @param transitionMat 
#' @param rateFunc 
#' @param tf 
#' @param replicates 
#' @param raw 
#' @param Ntf 
#' @param stats 
#'
#' @return a \code{Sim.ABM} object.
#' @examples Sim.ABM()
#' 
#' @export
setGeneric("Sim.ABM", function(params, N0, transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=15, 
                               raw=TRUE, discretePopSim=TRUE, Ntf=TRUE, stats=TRUE, ...) standardGeneric("Sim.ABM"))
setMethod("Sim.ABM",
          signature(params="list", N0="ANY", transitionsFunc="ANY", tf="ANY", replicates="ANY", raw="ANY", discretePopSim="ANY", Ntf="ANY", stats="ANY"),
          function(params){
            sim<- new("Sim.ABM", params=params)
            
            return (sim)
          }
)

setMethod("Sim.ABM",
          signature(params="missing", N0="ANY", transitionsFunc="ANY", tf="ANY", replicates="ANY", raw="ANY", discretePopSim="ANY", Ntf="ANY", stats="ANY"),
          function(N0, transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=15, raw=TRUE, discretePopSim=TRUE, Ntf=TRUE, stats=TRUE, ...){
            if (missing(N0)){
              N0<- c(N1s=0, N1b=1, N1bF=0, N2s=0, N2b=1, N2bF=0)
              N0<- lapply(2^(0:5), function(x) N0 * x)
              names(N0)<- paste0("N", sapply(N0, sum))
            }
            
            params<- list(N0=N0, transitionsFunc=transitionsFunc, tf=tf, replicates=replicates, raw=raw, discretePopSim=discretePopSim, Ntf=Ntf, stats=stats)
            
            dots<- list(...)
            params<- c(params, dots)
            
            sim<- Sim.ABM(params=params)
            
            return (sim)
          }
)


## Sim.ssa Class ----

#' @rdname Sim.ssa
#'
#' @param N0 
#' @param transitionMat 
#' @param rateFunc 
#' @param tf 
#' @param replicates 
#' @param raw 
#' @param Ntf 
#' @param stats 
#'
#' @return a \code{Sim.ssa} object.
#' @examples Sim.ssa()
#' 
#' @export
setGeneric("Sim.ssa", function(params, N0, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, 
                               tf=10, replicates=15, raw=TRUE, Ntf=TRUE, stats=TRUE, ...) standardGeneric("Sim.ssa"))
setMethod("Sim.ssa",
          signature(params="missing", N0="ANY", transitionMat="ANY", rateFunc="ANY", tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY", stats="ANY"),
          function(N0, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, tf=10, replicates=15, raw=TRUE, Ntf=TRUE, stats=TRUE, ...){
            if (missing(N0)){
              N0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
              N0<- lapply(2^(0:5), function(x) N0 * x)
              names(N0)<- paste0("N", sapply(N0, sum))
            }
            
            params<- list(N0=N0, transitionMat=transitionMat, rateFunc=rateFunc, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf, stats=stats)
            
            dots<- list(...)
            params<- c(params, dots)
            
            sim<- Sim.ssa(params=params)
            
            return (sim)
          }
)

setMethod("Sim.ssa",
          signature(params="list", N0="ANY", transitionMat="ANY", rateFunc="ANY", tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY", stats="ANY"),
          function(params){
            sim<- new("Sim.ssa", params=params)
            return (sim)
          }
)


## Generic ----
#' @export
setMethod("print", signature(x="Sim"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)

#' @export
setMethod("show", signature(object="Sim"),
          function(object){
            cat("Object of class", class(object), "with", nrow(object), "simulations\n Parameters:\n")
            cat("N0:\n")
            print(object@params$N0) # vector with more than one value
            tmp<- switch (class(object),
              Sim.ABM=print(data.frame(object@params[-c(1,2)]), row.names=FALSE), # one value only (1=N0, 2=transitionsFunc)
              Sim.ssa=print(data.frame(object@params[-c(1:3)]), row.names=FALSE), # one value only (1=N0, 2=transitionMat, 3=rateFunc)
              Sim=print(data.frame(object@params[-1]), row.names=FALSE) # one value only
            )
            
            cat("\n")
            print(S3Part(object))
            
          }
)


# Only allowed to subset by rows but $ and [[i]] works for columns
## Sim with results is ot useful without the complete Model. Use Model[] instead
# @rdname Sim
# @export
`[.Sim`<- function(x, ...){
  Sim(S3Part(x)[...])
}
