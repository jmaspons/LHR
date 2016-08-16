# Include functions for prototypes
# @include compoundDistributions.R

setOldClass("leslieMatrix")
setOldClass("discretePopSim")
setOldClass("numericDistri")
setOldClass("ssa")

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
setClass("Sim.discretePopSim", slots=list(Ntf="data.frame"), contains="Sim")
# TODO: , slots=list(raw="discretePopSim"), prototype=prototype(raw=discretePopSim())



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
# TODO: , slots=list(raw="numericDistri"), prototype=prototype(raw=distriBinom(1, .5))

#' SSA Simulation Class
#' 
#' @name Sim.ssa
#'
#' @slot Ntf data.frame. 
#'
#' @seealso \code{\link{Sim}}
#' @export
setClass("Sim.ssa", slots=list(Ntf="data.frame", deterministic="data.frame"), contains="Sim")
# TODO: , slots=list(raw="ssa"), prototype=prototype(raw=LHR:::exploreSSA())


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


#' Model.ssa
#' 
#' @name Model.ssa
#'
#' @export
setClass("Model.ssa", contains="Model")


## Check and overview ----
# getClass("ssa")
# getClass("Sim")

