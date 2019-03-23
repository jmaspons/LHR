numericDistri_default_options <- list(
  numericDistri.minP=1e-50
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(numericDistri_default_options) %in% names(op))
  if (any(toset)) options(numericDistri_default_options[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Expensive operations with numericDistri objects discard probabilities < ",
                        getOption("numericDistri.minP"),
                        ".\n Use: options(list(numericDistri.minP=1e-100) to change it.\n")
}