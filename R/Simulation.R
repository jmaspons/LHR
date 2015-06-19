## Simulation Class
# setClass("Sim", slots=list(discreteTime="discretePopSim"), contains="data.frame")
setClass("Sim", slots=list(params="list", raw="list"), 
         contains="data.frame")

## Constructor ----
setGeneric("Sim", function(params) standardGeneric("Sim"))

setMethod("Sim",
          signature(params="missing"),
          function(){
            Sim.discretePopSim()
          }
)

## Subclasses
## Sim.discretePopSim Class ----
setClass("Sim.discretePopSim", slots=list(Ntf="data.frame"), contains="Sim")
## TODO: don't export Sim.discretePopSim_complete
setGeneric("Sim.discretePopSim_complete", function(N0, envVar, sexRatio, matingSystem, tf, replicates, raw, Ntf) standardGeneric("Sim.discretePopSim_complete"))
setMethod("Sim.discretePopSim_complete",
          signature(N0="numeric", envVar="list", sexRatio="numeric", matingSystem="character", tf="numeric", replicates="numeric", raw="logical", Ntf="logical"),
          function(N0, envVar, sexRatio, matingSystem, tf, replicates, raw, Ntf){
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf)
            sim<- new("Sim.discretePopSim", params=params)
            return (sim)
          }
)

setGeneric("Sim.discretePopSim", function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=NA_real_, matingSystem=NA_character_, tf=10, replicates=15, raw=TRUE, Ntf=TRUE) standardGeneric("Sim.discretePopSim"))
setMethod("Sim.discretePopSim",
          signature(N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY", tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=NA_real_, matingSystem=NA_character_, tf=10, replicates=15, raw=TRUE, Ntf=TRUE){
            sim<- Sim.discretePopSim_complete(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf)
            return (sim)
          }
)


## Sim.numericDistri Class ----
#TODO: sexRatio
setClass("Sim.numericDistri", contains="Sim")
## TODO: don't export Sim.numericDistri_complete
setGeneric("Sim.numericDistri_complete", function(N0, envVar, sexRatio, matingSystem, tf, raw) standardGeneric("Sim.numericDistri_complete"))
setMethod("Sim.numericDistri_complete",
          signature(N0="numeric", envVar="list", sexRatio="numeric", matingSystem="character", tf="numeric", raw="logical"),
          function(N0, envVar, sexRatio, matingSystem, tf, raw){
            params<- list(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, raw=raw)
            sim<- new("Sim.numericDistri", params=params)
            return (sim)
          }
)

setGeneric("Sim.numericDistri", function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=NA_real_, matingSystem=NA_character_, tf=1, raw=TRUE) standardGeneric("Sim.numericDistri"))
setMethod("Sim.numericDistri",
          signature(N0="ANY", envVar="ANY", sexRatio="ANY", matingSystem="ANY", tf="ANY", raw="ANY"),
          function(N0=c(2, 10), envVar=list(j=TRUE, breedFail=FALSE), sexRatio=NA_real_, matingSystem=NA_character_, tf=1, raw=TRUE){
            sim<- Sim.numericDistri_complete(N0=N0, envVar=envVar, sexRatio=sexRatio, matingSystem=matingSystem, tf=tf, raw=raw)

            return (sim)
          }
)


## Sim.ssa Class ----
# setOldClass("ssa")
setClass("Sim.ssa", slots=list(Ntf="data.frame"), contains="Sim") #slots=list(result="ssa"), needs S3 class prototype
# TODO: don't export Sim.ssa_complete
setGeneric("Sim.ssa_complete", function(N0, transitionMat, rateFunc, tf, replicates, raw, Ntf, stats) standardGeneric("Sim.ssa_complete"))
setMethod("Sim.ssa_complete",
          signature(N0="list", transitionMat="closure", rateFunc="closure", tf="numeric", replicates="numeric", raw="logical", Ntf="logical", stats="logical"),
          function(N0, transitionMat, rateFunc, tf, replicates, raw, Ntf, stats){
            params<- list(N0=N0, transitionMat=transitionMat, rateFunc=rateFunc, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf, stats=stats)
            sim<- new("Sim.ssa", params=params)
            return (sim)
          }
)

setGeneric("Sim.ssa", function(N0, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, 
                               tf=10, replicates=15, raw=TRUE, Ntf=TRUE, stats=TRUE) standardGeneric("Sim.ssa"))
setMethod("Sim.ssa",
          signature(N0="ANY", transitionMat="ANY", rateFunc="ANY", tf="ANY", replicates="ANY", raw="ANY", Ntf="ANY", stats="ANY"),
          function(N0, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, tf=10, replicates=15, raw=TRUE, Ntf=TRUE, stats=TRUE){
            if (missing(N0)){
              N0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
              N0<- lapply(2^(0:5), function(x) N0 * x)
            }
            sim<- Sim.ssa_complete(N0=N0, transitionMat=transitionMat, rateFunc=rateFunc, tf=tf, replicates=replicates, raw=raw, Ntf=Ntf, stats=stats)
            return (sim)
          }
)



## Generic ----
setMethod("print", signature(x="Sim"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)

setMethod("show", signature(object="Sim"),
          function(object){
            cat("Object of class \"Sim\" with", nrow(object), "simulations\n Parameters:\n")
            print(object@params[1]) # vector with more than one value
            print(data.frame(object@params[-1]), row.names=FALSE) # one value only
            cat("\n")
            print(S3Part(object))
          }
)

setMethod("show", signature(object="Sim.ssa"),
          function(object){
            cat("Object of class \"Sim\" with", nrow(object), "simulations\n Parameters:\n")
            print(object@params[1]) # vector with more than one value
            print(data.frame(object@params[4:8]), row.names=FALSE) # one value only
            cat("\n")
            print(S3Part(object))
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.Sim`<- function(x, ...){
  Sim(S3Part(x)[...])
}