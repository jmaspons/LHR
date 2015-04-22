## Beta distribution parameter space ----
paramSpaceBeta<- function(method=c("meanVar", "freeShape"), pattern=c("montecarlo","regular"), n=30){
  pattern<- pattern[1]
  param<- switch(method[1], 
    meanVar={
      if (pattern == "regular"){
        betaDist<- seq(.05, .95, length=n)
        betaDist<- expand.grid(mean=betaDist, var=betaDist[which(betaDist<.25)])
      }else if (pattern == "montecarlo"){
        betaDist<- data.frame(mean=runif(n^2, .05, .95), var=runif(n^2,.05, .2499999999999999999999))
      }
      betaDist<- data.frame(betaDist, fbeta(betaDist$mean, betaDist$var))
      betaDist<- betaDist[!is.na(betaDist$shape1),]
    },
    freeShape={
      # runif
      ## 509 > a for Beta binomial
      ## 137 > a for Beta negative binomial
      if (pattern == "regular"){
        betaDist<- seq(1e-999, 100, length=n)
        betaDist<- expand.grid(shape1=betaDist, shape2=betaDist)
      }else if (pattern == "montecarlo"){
        betaDist<- data.frame(shape2=runif(n^2,1e-999, 100), shape1=runif(n^2,1e-999, 100))
      }
      betaDist<- data.frame(sbeta(betaDist$shape1, betaDist$shape2), betaDist)
    })

  return (param)
}

# For high dimensionality (e.g. interaction between 2 beta parameter space)
#   require(lhs)
#   hipercub<- maximinLHS(nP*1, 4)
#   betaDist2<- data.frame(shape1a=Escala(hipercub[,1], 1e-6, 100), shape2a=Escala(hipercub[,2], 1e-6, 100), shape1b=Escala(hipercub[,3], 1e-6, 100), shape2b=Escala(hipercub[,4], 1e-6, 100))

# tmp<- sbeta(betaDist$a, betaDist$b)