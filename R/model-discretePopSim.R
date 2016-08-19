#' Discrete Time Simulation
#' 
#' \code{discretePopSim} class represent the result of discrete time simulations with
#' replicates on rows and time on columns. The class inherits from \code{matrix} for 
#' subsetting by timesteps or replicates.
#' 
#' @name discretePopSim
#' @importFrom stats rbinom rbeta
#' @exportClass discretePopSim
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
setGeneric("discretePopSim", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, 
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
                    sexRatio="ANY", matingSystem="ANY", N0="ANY", replicates="ANY", tf="ANY", maxN="ANY"),
          function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0,
                   sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   N0=10, replicates=100, tf=10, maxN=100000){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }          
            
            broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE

            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0) varEnv<- TRUE       # environmental variation
            if (length(j) > 1 | length(a) > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes

            if (!broodMortality & !varEnv & !seasonalEnv & !sex){
              out<- mFit.t(fecundity=b, j=j, a=a, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & !seasonalEnv &  sex){
              out<- mFitSex.t(fecundity=b, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & seasonalEnv  & !sex){
              ## TODO: out<- mFit.tseason
            }
            if (!broodMortality & !varEnv & seasonalEnv  &  sex){
              ## TODO: out<- mFitSex.tseason
            }
            if (!broodMortality & varEnv  & !seasonalEnv & !sex){
              out<- mFit.tvar(fecundity=b, j=j, a=a, varJ=varJ, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & !seasonalEnv &  sex){
              out<- mFitSex.tvar(fecundity=b, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & seasonalEnv  & !sex){
              ## TODO: out<- mFit.tvarseason
            }
            if (!broodMortality & varEnv  & seasonalEnv  &  sex){
              ## TODO: out<- mFitSex.tvarseason
            }
            
            if (broodMortality  & !varEnv & !seasonalEnv & !sex){
              out<- mSurvBV.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & !seasonalEnv &  sex){
              out<- mSurvBVSex.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  & !sex){
              out<- mSurvBV.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  &  sex){
              out<- mSurvBVSex.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv & !sex){
              out<- mSurvBV.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv &  sex){
              out<- mSurvBVSex.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  & !sex){
              out<- mSurvBV.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  &  sex){
              out<- mSurvBVSex.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
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
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
  
  popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
  popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.tvar<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
mSurvBV.tseason<- function(broods, b, j, a, breedFail, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  for (t in 1:tf){
    for (k in seq_along(broods)){
      succeedingBroods<- rbinom(replicates, pop[,t], 1 - breedFail[k])
      # Juvenile survivors
      pop[,t+1]<- pop[,t+1] + rbinom(1, b * succeedingBroods, j[k])
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

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
## TODO
mSurvBV.tvarseason<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  ## TODO
  warning("Model not implemented?: mSurvBv.tvarseason")
  
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  if (varJ > 0) betaPars<- fbeta(mean=mean(j), var=varJ)
  if (varBreedFail > 0) breedFailEnv<- fbeta(mean=mean(breedFail), var=varBreedFail)
  jEnv<- j
  breedFailEnv<- breedFail
  
  for (t in 1:tf){
    # Environmental stochasticity
    if (varJ > 0){
      jEnv<- rbeta(replicates, shape1=betaPars$shape1, shape2=betaPars$shape2)
    }
    if (varBreedFail > 0){
      breedFailEnv<- rbeta(replicates, shape1=breedFailEnv$shape1, shape2=breedFailEnv$shape2)
    }
    
    #Seasonality no var
    if (length(j) > 1 & varJ == 0){
      
      succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - breedFailEnv)
      for (j in seq_along(succeedingBroods)){
        # Juvenile survivors
        pop[j,t+1]<- pop[j,t+1] + rbinom(1, b * succeedingBroods[j], sort(jEnv[1:succeedingBroods[j]], decreasing=TRUE))
      }
      
      #       envRep<- jEnv - mean(j) # Annual environment for each replicate
      #       jSeasonal<- merge(envRep, j)
      #       names(jSeasonal)<- c("envYear", "meanJseason")
      #       jSeasonal$j<- rowSums(jSeasonal)
      #       jSeasonal$replicate<- rep(1:replicates, times=length(j))
      
      
    }else{
      succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - breedFailEnv)
      # Juvenile survivors
      pop[,t+1]<- rbinom(replicates, b * succeedingBroods, jEnv)
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

  popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  
  betaPars<- fbeta(mean=j, var=varJ)
  
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    jEnv<- rbeta(replicates, shape1=jEnv$shape1, shape2=jEnv$shape2)
    
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
mSurvBVSex.tvar<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
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
mSurvBVSex.tseason<- function(broods, b, j, a, breedFail, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  
  for (t in 1:tf){
    for (k in seq_along(broods)){
      succeedingBroods<- rbinom(replicates, pop[,t], 1 - breedFail[k])
      # Juvenile survivors
      pop[,t+1]<- pop[,t+1] + rbinom(1, b * succeedingBroods, j[k])
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
mSurvBVSex.tvarseason<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  ## TODO
  warning("Model not implemented: mSurvBVSex.tvarseason")
  return (NA)
}
