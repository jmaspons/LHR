## Life History strategy class
# Contains sets of parameters defining the deterministic demographic parameters
# AFR: Age at First Reproduction
# a, j, s: survival for adults (A), juveniles [since egg/born until the first year] (J) and subadults [same as adults] (S)
# b: number of offspring per brood
# broods: number of broods per year
# fecundity = b  * broods: number of offspring per year
# fecundity * j: net fecundity

#' @include aaa-classes.R
NULL

## Constructors ----
#' @rdname LH
#'
#' @param pars 
#' @param lambda 
#' @param broods 
#' @param b clutch size
#' @param a adult survival
#' @param s subadult survival
#' @param j juvenile survival
#' @param AFR 
#' @param free which parameter is free to vary ("lambda", "j" or "a").
#' @param ... parameters passed to \code{\link{sampleLH}}.
#' 
#' @details Errors for free="a".
#' @return a \code{LH} object.
#' @examples 
#'  LH()
#'  LH(lambda=1)
#' 
#' @export
setGeneric("LH", function(pars, lambda=seq(.9, 1.1, by=0.1), fecundity, broods=2^(0:2), b=c(1, 2, 5, 10),
                          a=seq(0.3, 0.9, by=0.2), j=seq(0.2, 0.8, by=0.2), s=a, AFR=1, free="j", popbio=FALSE, ...) standardGeneric("LH"))

setMethod("LH",
          signature(pars="data.frame", lambda="missing", fecundity="missing", broods="missing", b="missing",
                    a="missing", j="missing", s="missing", AFR="missing", free="missing", popbio="ANY"),
          function(pars, popbio=FALSE){
            if (inherits(pars, "Model")) pars<- data.frame(pars, stringsAsFactors=TRUE)
            
            # if not defined, subadult survival is equal to adult survival. Only useful for AFR > 1
            if (!"s" %in% names(pars)){
              pars$s<- pars$a
            }
            
            if (!"idLH" %in% names(pars)) pars$idLH<- rownames(pars)
            
            selCols<- c("idLH", "lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR")
            if ("baseLH" %in% names(pars)) selCols<- c("baseLH", selCols)
            
            pars<- unique(pars[, selCols]) # Sort columns
            pars<- pars[naturalsort::naturalorder(pars$idLH),]
            rownames(pars)<- pars$idLH
            
            if (popbio & requireNamespace("popbio", quietly=TRUE)){
              popbio<- apply(pars[,-1], 1, function(x){ # First column is a character and makes x a character vector
                mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=fecundity * j, AFR=AFR))
                return(eigen.analisys2df(mat))
              })
              popbio<- do.call(rbind, popbio)
              
              pars<- cbind(pars, popbio)
            }
            
            return (new("LH", pars))
          }
)

setMethod("LH",
          signature(pars="missing", lambda="ANY", fecundity="missing", broods="ANY", b="ANY",
                    a="ANY", j="ANY", s="ANY", AFR="ANY", free="ANY", popbio="ANY"),
          function(lambda=seq(.9, 1.1, by=0.1), broods=2^(0:2), b=c(1, 2, 5, 10), 
                   a=seq(0.3, 0.9, by=0.2), j=seq(0.2, 0.8, by=0.2), s=a, AFR=1, free="j", popbio=FALSE, ...){

            pars<- sampleLH(lambda=lambda, broods=broods, b=b, j=j, a=a, AFR=AFR, free=free, ...)
            
            pars<- data.frame(idLH=rownames(pars), pars, stringsAsFactors=FALSE)
            
            LH(pars=pars, popbio=popbio)
          }
)


## Generic ----
#' @export
setMethod("print", signature(x="LH"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)

#' @export
setMethod("show", signature(object="LH"),
          function(object){
            cat("Object of class \"LH\" with", nrow(object), "strategies\n")
            print(S3Part(object))
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
#' @rdname LH
#' @export
`[.LH`<- function(x, ...){
            LH(pars=data.frame(x)[...])
}

## Sample LH ----

examplesLH<- function(){
  idLH<- c("fast", "slow", "freqRepro")
  lambda<- c(1.2, 1.05, 1.1)
  broods<- c(1, 1, 4)
  b<- c(10, 1, 1)
  a<- c(.4, .85, .6)
  s<- a
  AFR<- c(1, 1, 1) # ABM LH_behavior doesn't implement Age at First Reproduction
  
  fecundity<- broods * b
  
  j<- findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, AFR=AFR)
  
  pars<- data.frame(idLH, baseLH=idLH, lambda, fecundity, broods, b, a, s, j, AFR, stringsAsFactors=FALSE, row.names=idLH)
  
  return(pars)
}

# Impose the deterministic relations between the lambda and the rest of parameters
# lambda=1; broods=2^(0:2); b=1:10; j=c(.2, .8); a=seq(.3, .9, length=10); AFR=1
# sampleLH<- function(lambda=seq(.8, 2, by=0.1), broods=2^(0:2), b=c(1, seq(2, 20, by=2)), 
#                     j=seq(0.2, 0.8, by=0.1), a=seq(0.3, 0.9, by=0.1), AFR=1,
#                     free=c("j", "lambda")[1], maxFecundity=20, higherJuvMortality=TRUE, method=c("regular", "MonteCarlo"), census="pre-breeding"){

# WARNING("There are errors when estimating a = f(lambda, fecundity, j, AFR)")


#' Sample Life History strategies
#'
#' Impose the deterministic relations between the lambda and the rest of parameters
#' 
#' @param lambda 
#' @param broods 
#' @param b 
#' @param free 
#' @param maxFecundity 
#' @param higherJuvMortality 
#' @param method 
#' @param census 
#'
#' @return
#' @export
#'
#' @examples
sampleLH<- function(lambda=seq(1, 1.2, by=0.1), broods=2^(0:2), b=c(1, 2, 5, 10), 
                    j=seq(0.2, 0.8, by=0.2), a=seq(0.3, 0.9, by=0.2), AFR=1,
                    free=c("j", "lambda", "a"), maxFecundity=20, higherJuvMortality=TRUE, method=c("LH axes", "regular", "MonteCarlo"), census="pre-breeding"){
  free<- match.arg(free)
  method<- match.arg(method)
  
  if (method == "LH axes"){
    pars<- examplesLH()
    
    if (!missing(lambda)){
      comb<- expand.grid(lambda=lambda, idLH=pars$idLH)
      pars<- merge(pars[,-grep("lambda", names(pars))], comb, by="idLH")
      
      # Euler-Lotka corresponds to a pre-breding census matrix
      pars$j<- with(pars, findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, AFR=AFR))
      
      rownames(pars)<- paste0(pars$idLH, "-L", pars$lambda)
    }
    
    return(pars)
  }
  
  if (free == "lambda"){
    pars<- expand.grid(broods=broods, b=b, j=j, a=a, AFR=AFR)
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    
    # mean Lambda in the discrete time simulations correspons to pre-breeding census
    if (census == "pre-breeding"){
      pars$lambda<- apply(pars, 1, function(x){
        mat<- with(as.list(x), LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
        return(lambda(mat))
      })
    }else if (census == "post-breeding"){
      pars$lambda<- apply(pars, 1, function(x){
        mat<- with(as.list(x), LefkovitchPost(a=a, s=a, j=j, b=fecundity, AFR=AFR))  # subadult survival equal to adult survival # subadult survival equal to adult survival
        return(lambda(mat))
      })
    }
    
  }else if (free == "j"){
    pars<- expand.grid(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR)
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$j<- with(pars, findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, AFR=AFR))
  }else if (free == "a"){
    # warning("There are errors when estimating a = f(lambda, fecundity, j, AFR)")
    pars<- expand.grid(lambda=lambda, broods=broods, b=b, j=j, AFR=AFR)
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$a<- with(pars, findA_EulerLotka(lambda=lambda, b=fecundity, j=j, AFR=AFR))
  }
  
  pars<- stats::na.omit(pars)
  
  ## TODO: move to tests
  # Detect errors on the inverse eigenvalue problem and discard them
  if (free != "lambda"){
    lambdaMat<- apply(pars, 1, function(x){
      x<- as.list(x)
      if (is.na(x$j) | is.na(x$a)) return(NA)
      mat<- with(x, LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
      return(lambda(mat))
    })
    
    errLambda<- abs(pars$lambda - lambdaMat)
    errLambda<- which(errLambda > 0.001)
    
    if (length(errLambda) > 0) warning("Some errors on Euler-Lotka function to find parameters with ", free, "free.")
    
    if (free == "j"){
      pars$j[errLambda]<- NA
      pars<- pars[!is.na(pars$j),]
      pars<- pars[pars$j >= min(j),]
    }
    if (free == "a"){
      pars$a[errLambda]<- NA
      pars<- pars[!is.na(pars$a),]
      pars<- pars[pars$a >= min(a),]
    }
  }
  
  
  # Filter restrictions
  if (higherJuvMortality) pars<- pars[pars$j <= pars$a,]
  
  # Sort columns
  pars<- pars[order(pars$lambda, pars$a, pars$fecundity), c("lambda", "fecundity", "broods", "b", "a", "j", "AFR")]
  rownames(pars)<- NULL
  pars<- cbind(idLH=rownames(pars), pars)
  
  return (pars)
}
