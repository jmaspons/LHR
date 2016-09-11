#' Discrete Time Simulation
#' 
#' \code{discretePopSim} class represent the result of discrete time simulations with
#' replicates on rows and time on columns. The class inherits from \code{matrix} for 
#' subsetting by timesteps or replicates.
#' 
#' @name discretePopSim
#' @importFrom stats rbinom rbeta
NULL

# When the population gets extinct it fills results with NAs.
extinctNA<- function(pop){
  extT<- apply(pop, 1, function(x) match(0, x))
  extT[extT == ncol(pop)]<- NA 
  extPop<- which(!is.na(extT))
  for (i in extPop){
    pop[i,(extT[i]+1):ncol(pop)]<- NA
  }
  return(pop)
}

## S4 model dispatcher ----

#' Discrete time and population models
#'
#' TODO: add varA
#'
#' @rdname discretePopSim
#' @param broods 
#' @param b 
#' @param j 
#' @param a 
#' @param breedFail 
#' @param varJ 
#' @param varBreedFail 
#' @param seasonVar 
#' @param sexRatio 
#' @param matingSystem 
#' @param N0 
#' @param replicates 
#' @param tf 
#' @param maxN 
#'
#' @return a \code{discretePopSim} object.
#' @examples discretePopSim(b=2, j=.5, a=.5)
#'
#' @export
setGeneric("discretePopSim", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, seasonVar=1,
                                      sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), 
                                      N0=10, replicates=100, tf=10, maxN=100000) standardGeneric("discretePopSim"))
#### Common parameters:
# j=j, a=a, N0=N0, replicates=replicates, tf=tf, maxN=maxN
# 
#### Parameters by path:
# if(any(breedFail != 0)) # brood mortality: 
#   b=b, broods=broods, breedFail=breedFail
# else
#   b=b * broods
# 
# if (!is.na(matingSystem)) # two sexes: 
#   sexRatio=sexRatio, matingSystem=matingSystem
# 
# if (varJ != 0 | varBreedFail !=0) # environmental variation: 
#   varJ=varJ, varBreedFail=varBreedFail
# 
# if (length(j) > 1 | length(a) > 1) # environmental seasonality:
  
setMethod("discretePopSim", 
          signature(broods="ANY", b="ANY", j="numeric", a="numeric", breedFail="ANY", varJ="ANY", varBreedFail="ANY",
                    seasonVar="ANY", sexRatio="ANY", matingSystem="ANY", N0="ANY", replicates="ANY", tf="ANY", maxN="ANY"),
          function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, seasonVar=1,
                   sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   N0=10, replicates=100, tf=10, maxN=100000){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }          
            
            broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE
  
            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0) varEnv<- TRUE       # environmental variation
            if (any(seasonVar !=1) & broods > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes

            out<- NA

            if (!broodMortality & !varEnv & !seasonalEnv & !sex){
              out<- mFit.t(fecundity=b * broods, j=j, a=a, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & !seasonalEnv &  sex){
              out<- mFitSex.t(fecundity=b * broods, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & seasonalEnv  & !sex){
              out<- mFit.tseason(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
              ## NA: if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)
              # warning("## Shouldn't be called: out<- mFit.tseason\n",
              #         "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            if (!broodMortality & !varEnv & seasonalEnv  &  sex){
              ## NA: if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)
              warning("## Shouldn't be called: out<- mFitSex.tseason\n",
                      "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            if (!broodMortality & varEnv  & !seasonalEnv & !sex){
              out<- mFit.tvar(fecundity=b * broods, j=j, a=a, varJ=varJ, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & !seasonalEnv &  sex){
              out<- mFitSex.tvar(fecundity=b * broods, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & seasonalEnv  & !sex){
              out<- mFit.tvarseason(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, varJ=varJ, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
              ## NA: if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)
              # warning("## Shouldn't be called: out<- mFit.tvarseason\n",
              #         "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            if (!broodMortality & varEnv  & seasonalEnv  &  sex){
              ## NA: if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)
              warning("## Shouldn't be called: out<- mFitSex.tvarseason\n",
                      "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            
            if (broodMortality  & !varEnv & !seasonalEnv & !sex){
              out<- mSurvBV.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & !seasonalEnv &  sex){
              out<- mSurvBVSex.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  & !sex){
              out<- mSurvBV.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  &  sex){
              out<- mSurvBVSex.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv & !sex){
              out<- mSurvBV.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv &  sex){
              out<- mSurvBVSex.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  & !sex){
              out<- mSurvBV.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  &  sex){
              out<- mSurvBVSex.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            
            return (out)
          }
)


##  CONSTANT AND STABLE ENVIRONMENT ----
# Adult mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t, p=a) * fecundity, p=j)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
mFit.t<- function(fecundity, j, a, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, j)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.t<- function(broods, b, j, a, breedFail, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - breedFail)
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, b * succeedingBroods, j)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}


# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * fecundity, p=j), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
## TODO:
mFitSex.t<- function(fecundity, j, a, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  ## TODO
  warning("Model not implemented?: mFitSex.t")
  
  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    # Juvenile survivors
    reclutes<- rbinom(replicates, nPairs * fecundity * 2, j) # sex models should contemplate males (fecundity * 2)
    # Sex ratio
    popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
    popM[,t+1]<- reclutes - popF[,t+1]
    # Add adult survivors
    popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
    popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
    popF[which(popF[,t+1] > maxN),t+1]<- maxN
    popM[which(popM[,t+1] > maxN),t+1]<- maxN
  }
  
  popF<- popF[order(popF[,ncol(popF)]),]
  popM<- popM[order(popM[,ncol(popM)]),]
  popF<- extinctNA(popF)
  popM<- extinctNA(popM)
  class(popF)<- class(popM)<- c("discretePopSim", "matrix")
  
  pops<- list(females=popF, males=popM)
  
  return(pops)
}

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.t<- function(broods, b, j, a, breedFail, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    succeedingBroods<- rbinom(replicates, nPairs * broods, 1 - breedFail)
    # Juvenile survivors
    reclutes<- rbinom(replicates, succeedingBroods * b * 2, j)
    # Sex ratio
    popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
    popM[,t+1]<- reclutes - popF[,t+1]
    # Add adult survivors
    popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
    popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
    popF[which(popF[,t+1] > maxN),t+1]<- maxN
    popM[which(popM[,t+1] > maxN),t+1]<- maxN
  }
  
  popF<- popF[order(popF[,ncol(popF)]),]
  popM<- popM[order(popM[,ncol(popM)]),]
  popF<- extinctNA(popF)
  popM<- extinctNA(popM)
  class(popF)<- class(popM)<- c("discretePopSim", "matrix")
  
  pops<- list(females=popF, males=popM)
  
  return(pops)
}


## VARIABLE ENVIRONMENT ----
# j is a vector containing the different values for each brood attempt according to the seasonality (only for mSurvBV models)
## TODO: add seasonality and var combined
# survival<- data.frame(shape1=c(brood1, brood2...), shape2=c(brood1, brood2...))

# Adult mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t, p=a) * fecundity, p=Beta(j, varJ))
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
mFit.tvar<- function(fecundity, j, a, varJ, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  betaPars<- fbeta(mean=j, var=varJ)
  
  for (t in 1:tf){
    jEnv<- rbeta(replicates, shape1=betaPars$shape1, shape2=betaPars$shape2)
    
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, jEnv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

mFit.tseason<- function(broods, b, j, a, seasonVar, N0, replicates, tf, maxN=100000){
  jindSeason<- j * seasonVar
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  for (t in 1:tf){
    pop[,t+1]<- 0
    for (k in 1:broods){
      # Juvenile survivors
      pop[,t+1]<- pop[, t+1] + rbinom(replicates, b , jindSeason[k])
    }
    
    # Add adult survivors
    pop[, t+1]<- pop[, t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[, t+1] > maxN), t+1]<- maxN
  }
  
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

mFit.tvarseason<- function(broods, b, j, a, seasonVar, varJ=0, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  jindEnv<- j

  if (varJ > 0) betaJind<- fbeta(mean=j, var=varJ)

  for (t in 1:tf){
    jindEnv<- rbeta(replicates, shape1=betaJind$shape1, shape2=betaJind$shape2)
    
    pop[,t+1]<- 0
    for (k in 1:broods){
      # Juvenile survivors
      pop[,t+1]<- pop[, t+1] + rbinom(replicates, b, jindEnv * seasonVar[k])
    }
    
    # Add adult survivors
    pop[, t+1]<- pop[, t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[, t+1] > maxN), t+1]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.tvar<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  if (varJ > 0) betaParsJ<- fbeta(mean=j, var=varJ)
  if (varBreedFail > 0) betaParsBreedFail<- fbeta(mean=1 - breedFail, var=varBreedFail)
  jEnv<- j
  breedFailEnv<- breedFail
  
  for (t in 1:tf){
    # Environmental stochasticity
    if (varJ > 0){
      jEnv<- rbeta(replicates, shape1=betaParsJ$shape1, shape2=betaParsJ$shape2)
    }
    if (varBreedFail > 0){
      breedFailEnv<- rbeta(replicates, shape1=betaParsBreedFail$shape1, shape2=betaParsBreedFail$shape2)
    }
    
    succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - breedFailEnv)
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, b * succeedingBroods, jEnv)
    
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.tseason<- function(broods, b, j, a, breedFail, seasonVar, N0, replicates, tf, maxN=100000){
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[, 1]<- N0
  
  for (t in 1:tf){
    pop[,t+1]<- 0
    for (k in 1:broods){
      succeedingBroods<- rbinom(replicates, pop[, t], jbrSeason[k])
      # Juvenile survivors
      pop[,t+1]<- pop[, t+1] + rbinom(replicates, b * succeedingBroods, jindSeason[k])
    }
    
    # Add adult survivors
    pop[, t+1]<- pop[, t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[, t+1] > maxN), t+1]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
## TODO
mSurvBV.tvarseason<- function(broods, b, j, a, breedFail, seasonVar, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  jindEnv<- j
  jbrEnv<- 1 - breedFail
  
  if (varJ > 0) betaJind<- fbeta(mean=j, var=varJ)
  if (varBreedFail > 0) betaJbr<- fbeta(mean=breedFail, var=varBreedFail)
  
  for (t in 1:tf){
    # Environmental stochasticity
    if (varJ > 0){
      jindEnv<- rbeta(replicates, shape1=betaJind$shape1, shape2=betaJind$shape2)
    }
    
    if (varBreedFail > 0){
      jbrEnv<- rbeta(replicates, shape1=betaJbr$shape1, shape2=betaJbr$shape2)
    }
    
    pop[,t+1]<- 0
    
    for (k in 1:broods){
      succeedingBroods<- rbinom(replicates, pop[, t], jbrEnv * seasonVar[k])
      
      # Juvenile survivors
      pop[,t+1]<- pop[, t+1] + rbinom(replicates, b * succeedingBroods, jindEnv * seasonVar[k])
    }
      
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * fecundity, p=Beta(j, varJ)), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
## TODO:
mFitSex.tvar<- function(fecundity, j, a, varJ=0, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  ## TODO
  warning("Model not implemented?: mFitSex.tvar")

  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  
  betaPars<- fbeta(mean=j, var=varJ)
  
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    jEnv<- rbeta(replicates, shape1=betaPars$shape1, shape2=betaPars$shape2)
    
    # Juvenile survivors
    reclutes<- rbinom(replicates, nPairs * fecundity * 2, jEnv)
    # Sex ratio
    popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
    popM[,t+1]<- reclutes - popF[,t+1]
    # Add adult survivors
    popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
    popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
    popF[which(popF[,t+1] > maxN),t+1]<- maxN
    popM[which(popM[,t+1] > maxN),t+1]<- maxN
  }
  
  popF<- popF[order(popF[,ncol(popF)]),]
  popM<- popM[order(popM[,ncol(popM)]),]
  popF<- extinctNA(popF)
  popM<- extinctNA(popM)
  class(popF)<- class(popM)<- c("discretePopSim", "matrix")
  
  pops<- list(females=popF, males=popM)
  
  return(pops)
}

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.tvar<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry"), N0, replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)
  
  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  if (varJ > 0) betaParsJ<- fbeta(mean=j, var=varJ)
  if (varBreedFail > 0) betaParsBreedFail<- fbeta(mean=breedFail, var=varBreedFail)
  
  jEnv<- j
  breedFailEnv<- breedFail
  
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    if (varJ > 0){
      jEnv<- rbeta(replicates, shape1=betaParsJ$shape1, shape2=betaParsJ$shape2)
    }
    if (varBreedFail > 0){
      breedFailEnv<- rbeta(replicates, shape1=betaParsBreedFail$shape1, shape2=betaParsBreedFail$shape2)
    }
    
    succeedingBroods<- rbinom(replicates, nPairs * broods, 1 - breedFailEnv)
    # Juvenile survivors
    reclutes<- rbinom(replicates, succeedingBroods * b * 2, jEnv)
    # Sex ratio
    popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
    popM[,t+1]<- reclutes - popF[,t+1]
    # Add adult survivors
    popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
    popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
    popF[which(popF[,t+1] > maxN),t+1]<- maxN
    popM[which(popM[,t+1] > maxN),t+1]<- maxN
  }
  
  popF<- popF[order(popF[,ncol(popF)]),]
  popM<- popM[order(popM[,ncol(popM)]),]
  popF<- extinctNA(popF)
  popM<- extinctNA(popM)
  class(popF)<- class(popM)<- c("discretePopSim", "matrix")
  
  pops<- list(females=popF, males=popM)
  
  return(pops)
}


# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.tseason<- function(broods, b, j, a, breedFail, seasonVar, N0, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)

  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  pop<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  for (t in 1:tf){
    for (k in seq_along(broods)){
      succeedingBroods<- rbinom(replicates, pop[,t], jbrSeason[k])
      # Juvenile survivors
      pop[,t+1]<- pop[,t+1] + rbinom(1, b * succeedingBroods, jindSeason[k])
    }
    
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  
  pop<- pop[order(pop[,ncol(pop)]),]
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ)), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.tvarseason<- function(broods, b, j, a, breedFail, seasonVar, varJ=0, varBreedFail=0, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), N0, replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)
  
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  ## TODO
  warning("Model not implemented: mSurvBVSex.tvarseason")
  return (NA)
}

