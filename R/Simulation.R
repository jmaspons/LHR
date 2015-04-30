## Simulation Class
# setClass("Sim", slots=list(discreteTime="discretePopSim"), contains="data.frame")
setClass("Sim", slots=list(discreteTime="discretePopSim", probDistri="numericDistri", ssa="list"), contains="data.frame")

## Constructor ----
setGeneric("Sim", function(discret, distri, ssa) standardGeneric("Sim"))

setMethod("Sim",
          signature(discret="missing", distri="missing", ssa="missing"),
          function(ssa){
            simulation<- new("Sim")
            return (simulation)
          }
)

setMethod("Sim",
          signature(discret="discretePopSim", distri="missing", ssa="missing"),
          function(discret){
            simulation<- new("Sim", 
                            discret=discret)
            return (simulation)
          }
)

setMethod("Sim",
          signature(discret="missing", distri="numericDistri", ssa="missing"),
          function(distri){
            simulation<- new("Sim", 
                        distri=distri)
            return (simulation)
          }
)

setMethod("Sim",
          signature(discret="missing", distri="missing", ssa="ssa"),
          function(ssa){
            simulation<- new("Sim",
                        ssa=ssa)
            return (simulation)
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
            cat("Object of class \"Sim\" with", nrow(object), "simulations\n")
            print(S3Part(object))
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.Sim`<- function(x, ...){
  Sim(S3Part(x)[...])
}