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
## Sim.discretePopSim Class
setClass("Sim.discretePopSim", contains="Sim")
setGeneric("Sim.discretePopSim", function(N0, structure, tf, replicates, raw=TRUE) standardGeneric("Sim.discretePopSim"))
setMethod("Sim.discretePopSim",
          signature(N0="numeric", structure="character", tf="numeric", replicates="numeric", raw="logical"),
          function(N0, structure, tf, replicates, raw){
            params<- list(N0=N0, structure=structure, tf=tf, replicates=replicates, raw=raw)
            sim<- new("Sim.discretePopSim", 
                             params=params)
            return (sim)
          }
)

setMethod("Sim.discretePopSim",
          signature(N0="missing", structure="missing", tf="missing", replicates="missing", raw="missing"),
          function(){
            params<- list(N0=c(2, 10), structure="fit", tf=10, replicates=15, raw=TRUE)
            sim<- new("Sim.discretePopSim", 
                             params=params)
            return (sim)
          }
)

## Sim.numericDistri Class
setClass("Sim.numericDistri", contains="Sim")
setGeneric("Sim.numericDistri", function(N0, structure, raw) standardGeneric("Sim.numericDistri"))
setMethod("Sim.numericDistri",
          signature(N0="numeric", structure="character", raw="logical"),
          function(N0, structure, raw){
            params<- list(N0=N0, structure=structure, raw=raw)
            sim<- new("Sim.numericDistri", params=params)
            return (sim)
          }
)

setMethod("Sim.numericDistri",
          signature(N0="missing", structure="missing", raw="missing"),
          function(){
            params<- list(N0=c(2, 10), structure="fit", raw=TRUE)
            sim<- new("Sim.numericDistri", params=params)
            return (sim)
          }
)

## Sim.ssa Class
setClass("Sim.ssa", contains="Sim")
setGeneric("Sim.ssa", function(N0, structure, tf, replicates, raw) standardGeneric("Sim.ssa"))
setMethod("Sim.ssa",
          signature(N0="numeric", structure="character", tf="numeric", replicates="numeric", raw="logical"),
          function(N0, structure, tf, replicates, raw){
            params<- list(N0=N0, structure=structure, tf=tf, replicates, raw)
            simulation<- new("Sim.ssa", params=params)
            return (simulation)
          }
)

setMethod("Sim.ssa",
          signature(N0="missing", structure="missing", tf="missing", replicates="missing", raw="missing"),
          function(){
            params<- list(N0=c(2, 10), structure="LH-behavior", tf=10, replicates=15, raw=TRUE)
            sim<- new("Sim.ssa", params=params)
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

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.Sim`<- function(x, ...){
  Sim(S3Part(x)[...])
}