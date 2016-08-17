#' @include aaa-classes.R
NULL

## Constructors ----
#' @rdname LH.ssa
#' @param pars a data.frame.
#'
#' @export
setGeneric("LH.ssa", function(pars) standardGeneric("LH.ssa"))

setMethod("LH.ssa",
          signature(pars="data.frame"),
          function (pars){
            new("LH.ssa", pars)
          }
)
