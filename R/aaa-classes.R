
#' @exportClass discretePopSim
#' @exportClass discreteABMSim
#' @exportClass distriPopSim
#' @exportClass numericDistri
#' @exportClass leslieMatrix
NULL

setOldClass("discretePopSim")
setOldClass("discreteABMSim")
setOldClass("distriPopSim")
setOldClass("numericDistri")
setOldClass("leslieMatrix")


## Life history classes ----

#' Life histories class
#' 
#' @name LH
#'
#' @export
setClass("LH", contains="data.frame")


## Simulation classes ----

#' Simulation class
#' 
#' @name Sim
#'
#' @slot params list. 
#' @slot raw list. 
#' @slot N0_Pest.
#'
#' @export
setClass("Sim", slots=list(params="list", raw="ANY", N0_Pest="data.frame"), 
         contains="data.frame")


#' Discrete Simulation Class
#' 
#' @name Sim.discretePopSim
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.discretePopSim", slots=list(discretePopSim="list", Ntf="data.frame"), contains="Sim")
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

#' Model.discretePopSim
#' 
#' @name Model.discretePopSim
#'
#' @export
setClass("Model.discretePopSim", contains="Model")

#' Model.numericDistri
#' 
#' @name Model.numericDistri
#'
#' @export
setClass("Model.numericDistri", contains="Model")


#' Model.ABM
#' 
#' @name Model.ABM
#'
#' @export
setClass("Model.ABM", contains="Model")


## Check and overview ----
# getClass("Sim")

