## Life History strategy class
# Contains sets of parameters defining the deterministic demographic parameters
# AFR: Age at First Reproduction
# a, j, s: survival for adults (A), juveniles [since egg/born until the first year] (J) and subadults [same as adults] (S)
# b: number of offspring per brood
# broods: number of broods per year
# fecundity = b  * broods: number of offspring per year
# fecundity * j: net fecundity

setClass("LH", contains="data.frame")

## Constructors ----
setGeneric("LH", function(pars, lambda, fecundity, broods, b, a, s=a, j, AFR) standardGeneric("LH"))

setMethod("LH",
          signature(pars="data.frame", lambda="missing", fecundity="missing", broods="missing", b="missing",
                    a="missing", j="missing", s="missing", AFR="missing"),
          function(pars){
            # if not defined, subadult survival is equal to adult survival
            if (!"s" %in% names(pars)){
              pars$s<- pars$a
              pars<- pars[,c("lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR")]
            }
            strategy<- new("LH", pars)
            return (strategy)
          }
)

setMethod("LH",
          signature(pars="missing", lambda="missing", fecundity="missing", broods="missing", b="missing",
                    a="missing", j="missing", s="missing", AFR="missing"),
          function(){
            pars<- sampleLH()
            rownames(pars)<- NULL
            LH(pars)
          }
)

setMethod("LH",
          signature(pars="missing", lambda="numeric", fecundity="numeric", broods="numeric", b="numeric",
                    a="numeric", j="numeric", s="ANY", AFR="numeric"),
          function(lambda, fecundity, broods, b,   a, s=a, j, AFR=1){
            pars<- data.frame(lambda, fecundity, broods, b, a, s, j, AFR)
            strategy<- LH(pars)
            return (strategy)
          }
)

setMethod("LH",
          signature(pars="missing", lambda="missing", fecundity="numeric", broods="numeric", b="numeric",
                    a="numeric", j="numeric", s="ANY", AFR="numeric"),
          function(fecundity, broods, b,   a, s=a, j, AFR=1){
            pars<- sampleLH(broods=broods, b=b, j=j, a=a, AFR=AFR, free="lambda", 
                            maxFecundity=9999, higherJuvMortality=FALSE, census="pre-breeding")
            strategy<- LH(pars)
            return (strategy)
          }          
)

setMethod("LH",
          signature(pars="missing", lambda="numeric", fecundity="numeric", broods="numeric", b="numeric",
                    a="numeric", j="missing", s="ANY", AFR="numeric"),
          function(lambda, fecundity, broods, b,   a, s=a, AFR=1){
            pars<- sampleLH(lambda="numeric", broods=broods, b=b, a=a, AFR=AFR, free="j", 
                            maxFecundity=9999, higherJuvMortality=FALSE, census="pre-breeding")
            strategy<- LH(pars)
            return (strategy)
          }          
)

setMethod("LH",
          signature(pars="missing", lambda="numeric", fecundity="numeric", broods="numeric", b="numeric",
                    a="numeric", j="missing", s="ANY", AFR="numeric"),
          function(lambda, fecundity, broods, b,   a, s=a, AFR=1){
            pars<- sampleLH(lambda="numeric", broods=broods, b=b, j=j, AFR=AFR, free="a", 
                            maxFecundity=9999, higherJuvMortality=FALSE, census="pre-breeding")
            strategy<- LH(pars)
            return (strategy)
          }          
)

## Generic ----
setMethod("print", signature(x="LH"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)

setMethod("show", signature(object="LH"),
          function(object){
            cat("Object of class \"LH\" with", nrow(object), "strategies\n")
            print(S3Part(object))
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.LH`<- function(x, ...){
            LH(S3Part(x)[...])
}

## Sample LH imposing the deterministic relations between the lambda and the rest of parameters
sampleLH<- function(lambda=seq(.8, 2, by=0.1), broods=2^(0:2), b=c(1, seq(2, 20, by=1)), 
                    j=seq(0.2, 0.8, by=0.1), a=seq(0.3, 0.9, by=0.1), AFR=1,
                    free=c("j", "lambda")[1], maxFecundity=20, higherJuvMortality=TRUE, method=c("regular", "MonteCarlo"), census="pre-breeding"){
  if ("lambda" == free){
    pars<- expand.grid(broods, b, j, a, AFR)
    names(pars)<- c("broods", "b", "j", "a", "AFR")
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    
    pars$lambda<- NA
    for (i in 1:nrow(pars)){
      # mean Lambda in the discrete time simulations correspons to pre-breeding census
      if (census == "pre-breeding"){
        mat<- with(pars[i,], LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
      }else if (census == "post-breeding"){
        mat<- with(pars[i,], LefkovitchPost(a=a, s=a, j=j, b=fecundity, AFR=AFR))  # subadult survival equal to adult survival
      }
      pars$lambda[i]<- lambda(mat)
    }
  }else if (free == "j"){
    pars<- expand.grid(lambda, broods, b, a, AFR)
    names(pars)<- c("lambda", "broods", "b", "a", "AFR")
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$j<- with(pars, findJ_EulerLotka(lambda=lambda, b=b, a=a, AFR=AFR))
    # Filter
    #     pars<- pars[pars$j <= pars$a,]
  }else if (free == "a"){
    pars<- expand.grid(lambda, broods, b, j, AFR)
    names(pars)<- c("lambda", "broods", "b", "j", "AFR")
    pars$fecundity<- pars$broods * pars$b
    pars<- pars[pars$fecundity <= maxFecundity,]
    # Euler-Lotka corresponds to a pre-breding census matrix
    pars$a<- with(pars, findA_EulerLotka(lambda=lambda, b=b, j=j, AFR=AFR))
  }
  
  # Detect errors on the inverse eigenvalue problem and discard errors
  if (free != "lambda"){
    lambdaMat<- numeric(nrow(pars))
    for (i in 1:nrow(pars)){
      if (is.na(pars$j[i]) | is.na(pars$a[i])) next
      mat<- with(pars[i,], LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
      lambdaMat[i]<- lambda(mat)
    }
    errLambda<- abs(pars$lambda - lambdaMat)
    errLambda<- which(errLambda > 0.01)
    if (free == "j"){
      pars$j[errLambda]<- NA
      pars<- pars[!is.na(pars$j),]
    }
    if (free == "a"){
      pars$a[errLambda]<- NA
      pars<- pars[!is.na(pars$a),]
    }
  }
  
  # Filter restrictions
  if (higherJuvMortality) pars<- pars[pars$j <= pars$a,]
  
  # Sort columns
  pars<- pars[c("lambda", "fecundity", "broods", "b", "a", "j", "AFR")]
  return (pars)
}
