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
#' @param AFR Age at first reproduction
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
setGeneric("LH", function(pars, lambda=c(1.05, 1.2), fecundity, broods=2^(0:2), b=c(1, 2, 5, 10),
                          a=seq(0.3, 0.9, by=0.2), j=seq(0.2, 0.8, by=0.2), s=a, AFR=1, free="j", popbio=FALSE, ...) standardGeneric("LH"))

setMethod("LH",
          signature(pars="data.frame", lambda="missing", fecundity="missing", broods="missing", b="missing",
                    a="missing", j="missing", s="missing", AFR="missing", free="missing", popbio="ANY"),
          function(pars, popbio=FALSE){
            if (inherits(pars, "Model")) pars<- data.frame(pars, stringsAsFactors=FALSE)
            
            # if not defined, subadult survival is equal to adult survival. Effects for AFR > 1
            if (!"s" %in% colnames(pars)){
              pars$s<- pars$a
            }
            
            cols<- c("baseLH", "idLH", "lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR")
            colsPopbio<- c("elasFecundity", "elasSurvRepro", "elasSurvNonRepro", "generation.time", "net.reproductive.rate", "matureLifeExpectancy", "damping.ratio")
            
            if (popbio){
              cols<- c(cols, colsPopbio)
            }
            
            selCols<- intersect(cols, colnames(pars))
            pars<- unique(pars[, selCols])
            
            if (!"idLH" %in% colnames(pars)) pars<- data.frame(idLH=rownames(pars), pars, stringsAsFactors=FALSE)
            
            
            # Sort rows
            if ("baseLH" %in% names(pars)){
              pars<- pars[order(pars$baseLH, pars$lambda),]
            } else {
              pars<- pars[naturalsort::naturalorder(pars$idLH),]
            }
            
            rownames(pars)<- pars$idLH
            
            if (popbio){
              if (all(colsPopbio %in% colnames(pars))){
                reuse<- TRUE
              } else { reuse<- FALSE }
              
              if (!reuse & requireNamespace("popbio", quietly=TRUE)){
                # If one column is a character, makes x a character vector
                popbio<- apply(pars[, sapply(pars, is.numeric)], 1, function(x){
                  mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=fecundity * j, AFR=AFR))
                  return(eigen.analisys2df(mat))
                })
                
                popbio<- do.call(rbind, popbio)
                pars<- cbind(pars, popbio)
              }
            }
            
            return (new("LH", pars))
          }
)

setMethod("LH",
          signature(pars="missing", lambda="ANY", fecundity="missing", broods="ANY", b="ANY",
                    a="ANY", j="ANY", s="ANY", AFR="ANY", free="ANY", popbio="ANY"),
          function(lambda=c(1.05, 1.2), broods=2^(0:2), b=c(1, 2, 5, 10), 
                   a=seq(0.3, 0.9, by=0.2), j=seq(0.2, 0.8, by=0.2), s=a, AFR=1, free="j", popbio=FALSE, ...){

            pars<- sampleLH(lambda=lambda, broods=broods, b=b, j=j, a=a, s=s, AFR=AFR, free=free, ...)
            
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


#' Plot LH
#'
#' @rdname LH
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom graphics plot
plot.LH<- function(x, ...){
  x<- S3Part(x)
  if ("baseLH" %in% names(x)){
    x$colorLH<- factor(x$baseLH)
  } else {
    x$colorLH<- 1
  }
  
  cols<- intersect(names(x), c("lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR", "colorLH"))
  x<- unique(x[, cols])
  out<- graphics::plot(x[, sapply(x, is.numeric)], col=x$colorLH, ...) # All selected columns except colorLH
  # graphics::legend("topright", legend=levels(res$colorLH), bty = "y", pch = 19, col=res$colorLH)
  
  return(invisible(out))
}



## Sample LH ----

examplesLH<- function(){
  idLH<- c("fast", "slow", "freqRepro")
  lambda<- c(1.2, 1.05, 1.1)
  broods<- c(1, 1, 4)
  b<- c(10, 2, 1)
  a<- c(.4, .85, .6)
  s<- c(.4, .85, .6) # only used in AFR > 1
  AFR<- c(1, 4, 1)
  
  fecundity<- broods * b
  
  j<- findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, s=s, AFR=AFR)
  
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
                    j=seq(0.2, 0.8, by=0.2), a=seq(0.3, 0.9, by=0.2), s, AFR=2^(0:2),
                    free=c("j", "lambda", "a"), maxFecundity=20, higherJuvMortality=TRUE, method=c("LH axes", "regular", "MonteCarlo"), census="pre-breeding"){
  free<- match.arg(free)
  method<- match.arg(method)
  
  if (method == "LH axes"){
    pars<- examplesLH()
    
    if (!missing(lambda)){
      comb<- expand.grid(lambda=lambda, idLH=pars$idLH, stringsAsFactors=FALSE)
      pars<- merge(pars[,-grep("lambda", names(pars))], comb, by="idLH")
      
      # Euler-Lotka corresponds to a pre-breding census matrix
      pars$j<- with(pars, findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, s=s, AFR=AFR))
      
      pars$idLH<- paste0(pars$idLH, "-L", pars$lambda)
      rownames(pars)<- pars$idLH
      pars<- pars[order(pars$baseLH, pars$lambda),]
    }
    
    if (any(is.na(pars))){
      warning("Some parameter combinations produce NAs and are discarded (eg. probabilities > 1)")
      pars<- stats::na.omit(pars)
    }
    
    return(pars)
  }
  
  misS<- missing(s) # s == a

  if (free == "lambda"){
    if (misS){
      pars<- expand.grid(broods=broods, b=b, j=j, a=a, AFR=AFR)
      pars$s<- pars$a
    }else{
      pars<- expand.grid(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR)
    }
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    
    # mean Lambda in the discrete time simulations correspons to pre-breeding census
    if (census == "pre-breeding"){
      pars$lambda<- apply(pars, 1, function(x){
        mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=j * fecundity, AFR=AFR))
        return(lambda(mat))
      })
    }else if (census == "post-breeding"){
      pars$lambda<- apply(pars, 1, function(x){
        mat<- with(as.list(x), LefkovitchPost(a=a, s=s, j=j, b=fecundity, AFR=AFR))
        return(lambda(mat))
      })
    }
    
  }else if (free == "j"){
    if (misS){
      pars<- expand.grid(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR)
      pars$s<- pars$a
    }else{
      pars<- expand.grid(lambda=lambda, broods=broods, b=b, s=s, a=a, AFR=AFR)
    }
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$j<- with(pars, findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, s=s, AFR=AFR))
  }else if (free == "a"){
    if (misS){
      pars<- expand.grid(lambda=lambda, broods=broods, b=b, j=j, AFR=AFR)
      pars$s<- pars$j ## WARNING: s<- j
    }else{
      pars<- expand.grid(lambda=lambda, broods=broods, b=b, j=j, s=s, AFR=AFR)
    }
    # warning("There are errors when estimating a = f(lambda, fecundity, j, AFR)")
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$a<- with(pars, findA_EulerLotka(lambda=lambda, b=fecundity, j=j, s=s, AFR=AFR))
  }
  
  if (any(is.na(pars))){
    warning("Some parameter combinations produce NAs and are discarded (eg. probabilities > 1)")
    pars<- stats::na.omit(pars)
  }
  
  # Filter restrictions
  if (higherJuvMortality) pars<- pars[pars$j <= pars$a,]
  
  # Sort columns
  pars<- pars[order(pars$lambda, pars$a, pars$fecundity), c("lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR")]
  rownames(pars)<- NULL
  pars<- cbind(idLH=rownames(pars), pars, stringsAsFactors=FALSE)
  
  return (pars)
}
