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
setGeneric("Sim", function(params, type=c("discretePopSim", "numericDistri", "ABM", "ssa")) standardGeneric("Sim"))

setMethod("Sim",
          signature(params="missing", type="ANY"),
          function(params, type=c("discretePopSim", "numericDistri", "ABM", "ssa")){
            type<- match.arg(type)
            
            sim<- switch(type,
                         discretePopSim=Sim.discretePopSim(),
                         numericDistri=Sim.numericDistri(),
                         ABM=Sim.ABM(),
                         ssa=Sim.ssa())
            return (sim)
          }
)

setMethod("Sim",
          signature(params="list", type="character"),
          function(params, type=c("discretePopSim", "numericDistri", "ABM", "ssa")){
            type<- match.arg(type)
            
            sim<- switch(type,
                   discretePopSim=Sim.discretePopSim(params=params),
                   numericDistri=Sim.numericDistri(params=params),
                   ABM=Sim.ABM(params=params),
                   ssa=Sim.ssa(params=params))

            return (sim)
          }
)

setMethod("Sim",
          signature(params="Model", type="missing"),
          function(params){
            
            type<- class(params@sim)
            simPars<- params@sim@params
            
            sim<- switch(type,
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
setGeneric("Sim.discretePopSim", function(params, N0=c(2, 10), envVar=list(j=TRUE, breedFail=TRUE),
                                          sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                                          tf=10, replicates=15, maxN=10000, raw=TRUE, Ntf=TRUE) standardGeneric("Sim.discretePopSim"))

setMethod("Sim.discretePopSim",
          signature(params="list", N0="missing", envVar="missing", sexRatio="missing", matingSystem="missing", tf="missing", replicates="missing", maxN="missing", raw="missing", Ntf="missing"),
          function(params){
            sim<- new("Sim.discretePopSim", params=params)

            return (sim)
          }
)

setMethod("Sim.discretePopSim",
          signature(params="missing", N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY",
                    tf="ANY", replicates="ANY", maxN="ANY", raw="ANY", Ntf="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=TRUE), sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   tf=10, replicates=15, maxN=10000, raw=TRUE, Ntf=TRUE){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }
            
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, replicates=replicates, maxN=maxN, raw=raw, Ntf=Ntf)
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
setGeneric("Sim.numericDistri", function(params, N0=c(2, 10), envVar=list(j=TRUE, breedFail=TRUE), 
                                         sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), 
                                         tf=1, maxN=10000, raw=TRUE) standardGeneric("Sim.numericDistri"))

setMethod("Sim.numericDistri",
          signature(params="list", N0="missing", envVar="missing", sexRatio="missing", matingSystem="missing", tf="missing", maxN="missing", raw="missing"),
          function(params){
            sim<- new("Sim.numericDistri", params=params)
            
            return (sim)
          }
)

setMethod("Sim.numericDistri",
          signature(params="missing", N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY", tf="ANY", maxN="ANY", raw="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=TRUE), sexRatio=0.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), tf=1, maxN=10000, raw=TRUE){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }
            
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, maxN=maxN, raw=raw)
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
setGeneric("Sim.ABM", function(params, N0, transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=15, maxN=10000,
                               raw=TRUE, discretePopSim=TRUE, Ntf=TRUE, stats=TRUE, ...) standardGeneric("Sim.ABM"))
setMethod("Sim.ABM",
          signature(params="list", N0="missing", transitionsFunc="missing", tf="missing", replicates="missing", maxN="missing", raw="missing", discretePopSim="missing", Ntf="missing", stats="missing"),
          function(params){
            sim<- new("Sim.ABM", params=params)
            
            return (sim)
          }
)

setMethod("Sim.ABM",
          signature(params="missing", N0="ANY", transitionsFunc="ANY", tf="ANY", replicates="ANY", maxN="ANY", raw="ANY", discretePopSim="ANY", Ntf="ANY", stats="ANY"),
          function(N0, transitionsFunc=transitionABM.LH_Beh, tf=10, replicates=15, maxN=10000, raw=TRUE, discretePopSim=TRUE, Ntf=TRUE, stats=TRUE, ...){
            if (missing(N0)){
              N0<- c(N1s=0, N1b=1, N1bF=0, N2s=0, N2b=1, N2bF=0)
              N0<- lapply(2^(0:5), function(x) N0 * x)
              names(N0)<- paste0("N", sapply(N0, sum))
            }
            
            params<- list(N0=N0, transitionsFunc=transitionsFunc, tf=tf, replicates=replicates, maxN=maxN, raw=raw, discretePopSim=discretePopSim, Ntf=Ntf, stats=stats)
            
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
                               tf=10, replicates=15, maxN=10000, raw=FALSE, discretePop=TRUE, Ntf=TRUE, stats=TRUE, ...) standardGeneric("Sim.ssa"))
setMethod("Sim.ssa",
          signature(params="missing", N0="ANY", transitionMat="ANY", rateFunc="ANY", tf="ANY", replicates="ANY", maxN="ANY", raw="ANY", discretePop="ANY", Ntf="ANY", stats="ANY"),
          function(N0, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, tf=10, replicates=15, maxN=10000, raw=FALSE, discretePop=TRUE, Ntf=TRUE, stats=TRUE, ...){
            if (missing(N0)){
              N0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
              N0<- lapply(2^(0:5), function(x) N0 * x)
              names(N0)<- paste0("N", sapply(N0, sum))
            }
            
            params<- list(N0=N0, transitionMat=transitionMat, rateFunc=rateFunc, tf=tf, replicates=replicates, maxN=maxN, raw=raw, discretePop=discretePop, Ntf=Ntf, stats=stats)
            
            dots<- list(...)
            params<- c(params, dots)
            
            sim<- Sim.ssa(params=params)
            
            return (sim)
          }
)

setMethod("Sim.ssa",
          signature(params="list", N0="missing", transitionMat="missing", rateFunc="missing", tf="missing", replicates="missing", maxN="missing", raw="missing", discretePop="missing", Ntf="missing", stats="missing"),
          function(params){
            sim<- new("Sim.ssa", params=params)
            return (sim)
          }
)


## Generic ----

#' @export
setMethod("show", signature(object="Sim"),
          function(object){
            cat("Object of class", class(object), "with", nrow(object), "simulations\n\n")
            
            if (nrow(object) > 0 ){
              print(S3Part(object))
              cat("\n\n")
            }
            
            if (nrow(object@N0_Pest) > 0 ){
              cat("N0=f(P_est):\n")
              print(object@N0_Pest)
              cat("\n\n")
            }
            
            cat(" Parameters:\n")
            pars<- object@params
            
            func<- sapply(pars, inherits, "function")
            gt1<- sapply(pars, length) > 1
            
            print(pars[gt1 & !func]) # Parameters with more than one value
            print(data.frame(pars[!gt1 &!func]), row.names=FALSE) # Parameters with only one value
            
            invisible(object)
          }
)


# Only allowed to subset by rows but $ and [[i]] works for columns
## Sim with results is ot useful without the complete Model. Use Model[] instead
# @rdname Sim
# @export
`[.Sim`<- function(x, ...){
  Sim(S3Part(x)[...])
}
