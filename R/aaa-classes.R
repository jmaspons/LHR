# Include functions for prototypes
# @include compoundDistributions.R

#' @exportClass discretePopSim
#' @exportClass discreteABMSim
#' @exportClass exploreABMSim
#' @exportClass numericDistri
#' @exportClass ssa
#' @exportClass leslieMatrix
NULL

setOldClass("discretePopSim")
setOldClass("discreteABMSim")
setOldClass("exploreABMSim")
setOldClass("numericDistri")
setOldClass("ssa")
setOldClass("leslieMatrix")


## Life history classes ----

#' Life histories class
#' 
#' @name LH
#'
#' @export
setClass("LH", contains="data.frame")


#' LH.ssa
#' 
#' @name LH.ssa
#'
#' @export
setClass("LH.ssa", contains="LH")


## Simulation classes ----

# setClass("Sim", slots=list(discreteTime="discretePopSim"), contains="data.frame")
#' Simulation class
#' 
#' @name Sim
#'
#' @slot params list. 
#' @slot raw list. 
#'
#' @export
setClass("Sim", slots=list(params="list", raw="ANY"), 
         contains="data.frame")


#' Discrete Simulation Class
#' 
#' @name Sim.discretePopSim
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.discretePopSim", slots=list(discretePopSim="ANY", Ntf="data.frame"), contains="Sim")
# TODO: , slots=list(discretePopSim="discretePopSim"), prototype=prototype(discretePopSim=discretePopSim())



#TODO: sexRatio
#' Numeric Distribution Simulation Class
#' 
#' @name Sim.numericDistri
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.numericDistri", contains="Sim")

#' Agent Based Model Simulation Class
#' 
#' @name Sim.ABM
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.ABM", slots=list(deterministic="data.frame"), contains="Sim.discretePopSim")

#' SSA Simulation Class
#' 
#' @name Sim.ssa
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.ssa", slots=list(deterministic="data.frame"), contains="Sim.discretePopSim")


## Environment class ----

#' Environment class
#' 
#' @name Env
#' @slot seasonRange data.frame. 
#'
#' @examples Env()
#' @export
setClass("Env", slots=list(seasonRange="data.frame"), contains="data.frame")


## Model classes ----

#' Model class
#'
#' @name Model
#' 
#' @slot sim Sim. 
#' @slot params list. 
#'
#' @export
setClass("Model", slots=list(sim="Sim", params="list"), contains="data.frame")

#' Model.ABM
#' 
#' @name Model.ABM
#'
#' @export
setClass("Model.ABM", contains="Model")

#' Model.ssa
#' 
#' @name Model.ssa
#'
#' @export
setClass("Model.ssa", contains="Model")


## Check and overview ----
# getClass("ssa")
# getClass("Sim")

