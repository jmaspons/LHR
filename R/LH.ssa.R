#' LH.ssa
#'
#' @include LH.R
#' @export
setClass("LH.ssa", contains="LH")

## Constructors ----
#' @describeIn LH.ssa
#' @param pars a data.frame.
#'
#' @export
setGeneric("LH.ssa", function(pars) standardGeneric("LH.ssa"))

setMethod("LH.ssa",
          signature(pars="data.frame"),
          function (pars){
            new("LH.ssa")
          }
)
