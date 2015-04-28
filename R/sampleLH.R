## Parameters
# AFR: Age at First Reproduction
# a, j, s: survival for adults (A), juveniles [since egg/born until the first year] (J) and subadults [same as adults] (S)
# b: number of offspring per brood
# broods: number of broods per year
# fecundity = b  * broods: number of offspring per year
# fecundity * j: net fecundity

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

