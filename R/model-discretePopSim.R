## TODO: add subadult classes for models with 2 sexes and AFR > 1
## TODO: add varA for *Sex.tvar models

#' Discrete Time Simulation
#' 
#' \code{discretePopSim} class represent the result of discrete time simulations with
#' replicates on rows and time on columns. The class inherits from \code{matrix} for 
#' subsetting by timesteps or replicates.
#' 
#' @name discretePopSim
#' @importFrom stats rbinom rbeta
NULL

cleanDiscretePopSim<- function(pop, maxN){
  if (ncol(pop) == 2) return(pop) ## for ABM with Ntf=TRUE (keep N0 and Ntf only)
  
  pop<- extinctNA(pop)
  pop<- maxNNA(pop, maxN=maxN)

  return(pop)
}

# Fill results with NAs after the population gets extinct (pop == 0).
extinctNA<- function(pop){
  extT<- apply(pop, 1, function(x) match(0, x))
  extT[which(extT == ncol(pop))]<- NA 
  extPop<- which(!is.na(extT))
  
  tmp<- cbind(extT=extT[extPop], pop[extPop, , drop=FALSE])
  
  tmp<- apply(tmp, 1, function(x){
    x[(x[1]+1):ncol(pop) + 1]<- NA
    x[-1]
  })
  pop[extPop, ]<- t(tmp)
  
  return(pop)
}

# Fill results with NAs after the population reach maxN
maxNNA<- function(pop, maxN){
  if (missing(maxN)){
    maxN<- max(pop[, -1], na.rm=TRUE) # 1 : N0
    if (!any(apply(pop[, -1, drop=FALSE], 2, function(x) all(x == maxN)), na.rm=TRUE)){ # not all populations reach maxN
      return(pop)
    }
  }

  maxNT<- apply(pop[, -1, drop=FALSE], 1, function(x) match(maxN, x)) + 1
  maxNT[which(maxNT == ncol(pop))]<- NA 
  maxNPop<- which(!is.na(maxNT))
  
  tmp<- cbind(maxNT=maxNT[maxNPop], pop[maxNPop, , drop=FALSE])
  
  tmp<- apply(tmp, 1, function(x){
    x[(x[1]+1):ncol(pop) + 1]<- NA
    x[-1]
  })
  pop[maxNPop, ]<- t(tmp)
  
  return(pop)
}


## S4 model dispatcher ----

#' Discrete time and population models
#'
#' @rdname discretePopSim
#' @param broods number of broods per year
#' @param b offsprings per brood
#' @param j juvenile survival
#' @param s subadult survival
#' @param a adult survival
#' @param AFR age at first reproduction
#' @param breedFail 
#' @param varJ 
#' @param varBreedFail 
#' @param varA
#' @param seasonVar 
#' @param sexRatio 
#' @param matingSystem 
#' @param N0 
#' @param replicates number of populations to simulate
#' @param tf final time
#' @param maxN maximum population size
#'
#' @return a \code{discretePopSim} object.
#' @examples discretePopSim(b=2, j=.5, a=.5)
#'
#' @export
setGeneric("discretePopSim", function(broods=1, b, j, s=a, a, AFR=1, breedFail=0, varJ=0, varBreedFail=0, varA=0, seasonVar=1,
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
          signature(broods="ANY", b="ANY", j="numeric", s="ANY", a="numeric", AFR="ANY", breedFail="ANY", varJ="ANY", varBreedFail="ANY", varA="ANY",
                    seasonVar="ANY", sexRatio="ANY", matingSystem="ANY", N0="ANY", replicates="ANY", tf="ANY", maxN="ANY"),
          function(broods=1, b, j, s=a, a, AFR=1, breedFail=0, varJ=0, varBreedFail=0, varA=0, seasonVar=1,
                   sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   N0=10, replicates=100, tf=10, maxN=100000){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }          
            
            broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE
  
            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0 | varA !=0) varEnv<- TRUE       # environmental variation
            if (any(seasonVar !=1) & broods > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes

            out<- NA

            if (!broodMortality & !varEnv & !seasonalEnv & !sex){
              out<- mFit.t(fecundity=b * broods, j=j, s=s, a=a, AFR=AFR, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & !seasonalEnv &  sex){
              out<- mFitSex.t(fecundity=b * broods, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & seasonalEnv  & !sex){
              out<- mFit.tseason(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, seasonVar=seasonVar, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
              ## NA: if there is no broodMortality or only 1 brood, seasonality have no effect (always breed on th best period)
              # warning("## Shouldn't be called: out<- mFit.tseason\n",
              #         "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            if (!broodMortality & !varEnv & seasonalEnv  &  sex){
              ## NA: if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)
              warning("## Shouldn't be called: out<- mFitSex.tseason\n",
                      "if there is no broodMortality or only 1 brood seasonality have no effect (always breed on th best period)")
            }
            if (!broodMortality & varEnv  & !seasonalEnv & !sex){
              out<- mFit.tvar(fecundity=b * broods, j=j, s=s, a=a, AFR=AFR, varJ=varJ, varA=varA, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & !seasonalEnv &  sex){
              out<- mFitSex.tvar(fecundity=b * broods, j=j, a=a, varJ=varJ, varA=varA, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & seasonalEnv  & !sex){
              out<- mFit.tvarseason(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, seasonVar=seasonVar, varJ=varJ, varA=varA, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
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
              out<- mSurvBV.t(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & !seasonalEnv &  sex){
              out<- mSurvBVSex.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  & !sex){
              out<- mSurvBV.tseason(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, breedFail=breedFail, seasonVar=seasonVar, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  &  sex){
              out<- mSurvBVSex.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv & !sex){
              out<- mSurvBV.tvar(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, varA=varA, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv &  sex){
              out<- mSurvBVSex.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, varA=varA, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  & !sex){
              out<- mSurvBV.tvarseason(broods=broods, b=b, j=j, s=s, a=a, AFR=AFR, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, varA=varA, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  &  sex){
              out<- mSurvBVSex.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, varA=varA, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            
            return (out)
          }
)


##  CONSTANT AND STABLE ENVIRONMENT ----
# Adult mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t, p=a) * fecundity, p=j)
# fecundity = b x broods
mFit.t<- function(fecundity, j, s, a, AFR, N0, replicates, tf, maxN=100000){
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  for (t in 1:tf){
    # Juvenile survivors
    pop[, t+1, 1]<- rbinom(replicates, pop[, t, AFR] * fecundity, j)
    # Add adult survivors
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], a)
    # Subadults (growth transitions)
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, s))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, s))
    }
    
    pop[which(pop[,t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j)
mSurvBV.t<- function(broods, b, j, s, a, AFR, breedFail, N0, replicates, tf, maxN=100000){
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  for (t in 1:tf){
    succeedingBroods<- rbinom(replicates, pop[, t, AFR] * broods, 1 - breedFail)
    # Juvenile survivors
    pop[, t+1, 1]<- rbinom(replicates, b * succeedingBroods, j)
    # Add adult survivors
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], a)
    # Subadults (growth transitions)
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, s))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, s))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}


# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * fecundity, p=j), p=sexRatio)
# fecundity = b x broods
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# TODO: add subadults
mFitSex.t<- function(fecundity, j, s, a, AFR, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
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
# TODO: add subadults
mSurvBVSex.t<- function(broods, b, j, s, a, AFR, breedFail, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
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

# Adult mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t, p=a) * fecundity, p=Beta(j, varJ))
# fecundity = b x broods
mFit.tvar<- function(fecundity, j, s, a, AFR, varJ, varA, N0, replicates, tf, maxN=100000){
  
  if (varJ > 0){
    betaParsJ<- fbeta(mean=j, var=varJ)
    if (anyNA(betaParsJ)){
      warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  
  if (varA > 0){
    betaParsA<- fbeta(mean=a, var=varA)
    if (anyNA(betaParsA)){
      warning("The combination of adult survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
    
    if (AFR > 1){
      betaParsS<- fbeta(mean=s, var=varA)
      if (anyNA(betaParsS)){
        warning("The combination of subadult survival and variance fall outside the domain of the Beta distribution.")
        return (NA)
      }
    }
  }
  
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  jEnv<- j
  sEnv<- s
  aEnv<- a

  for (t in 1:tf){
    # Juvenile survivors
    if (varJ > 0)
      jEnv<- rbeta(replicates, shape1=betaParsJ$shape1, shape2=betaParsJ$shape2)

    pop[, t+1, 1]<- rbinom(replicates, pop[, t, AFR] * fecundity, jEnv)
    
    # Add adult survivors
    if (varA > 0)
      aEnv<- rbeta(replicates, shape1=betaParsA$shape1, shape2=betaParsA$shape2)
    
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], aEnv)
    
    # Subadults (growth transitions)
    if (varA > 0 & AFR > 1)
      sEnv<- rbeta(replicates, shape1=betaParsS$shape1, shape2=betaParsS$shape2)
    
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, sEnv))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, sEnv))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

mFit.tseason<- function(broods, b, j, s, a, AFR, seasonVar, N0, replicates, tf, maxN=100000){
  jindSeason<- j * seasonVar
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  for (t in 1:tf){
    for (k in 1:broods){
      # Juvenile survivors
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, pop[, t, AFR] * b, jindSeason[k])
    }
    
    # Add adult survivors
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], a)
    
    # Subadults (growth transitions)
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, s))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, s))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

mFit.tvarseason<- function(broods, b, j, s, a, AFR, seasonVar, varJ=0, varA=0, N0, replicates, tf, maxN=100000){
  
  if (varJ > 0){
    betaParsJ<- fbeta(mean=j, var=varJ)
    if (anyNA(betaParsJ)){
      warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }

  if (varA > 0){
    betaParsA<- fbeta(mean=a, var=varA)
    if (anyNA(betaParsA)){
      warning("The combination of adult survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
    
    if (AFR > 1){
      betaParsS<- fbeta(mean=s, var=varA)
      if (anyNA(betaParsS)){
        warning("The combination of subadult survival and variance fall outside the domain of the Beta distribution.")
        return (NA)
      }
    }
  }
  
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  jEnv<- j
  sEnv<- s
  aEnv<- a
  
  for (t in 1:tf){
    if (varJ > 0)
      jEnv<- rbeta(replicates, shape1=betaParsJ$shape1, shape2=betaParsJ$shape2)
    
    for (k in 1:broods){
      # Juvenile survivors
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, pop[, t, AFR] * b, jEnv * seasonVar[k])
    }
    
    # Add adult survivors
    if (varA > 0)
      aEnv<- rbeta(replicates, shape1=betaParsA$shape1, shape2=betaParsA$shape2)
    
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], aEnv)
    
    # Subadults (growth transitions)
    if (varA > 0 & AFR > 1)
      sEnv<- rbeta(replicates, shape1=betaParsS$shape1, shape2=betaParsS$shape2)
    
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, sEnv))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, sEnv))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
mSurvBV.tvar<- function(broods, b, j, s, a, AFR, breedFail, varJ=0, varBreedFail=0, varA=0, N0, replicates, tf, maxN=100000){
  
  if (varJ > 0){
    betaParsJ<- fbeta(mean=j, var=varJ)
    if (anyNA(betaParsJ)){
      warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  if (varBreedFail > 0){
    betaParsBreedFail<- fbeta(mean=breedFail, var=varBreedFail)
    if (anyNA(betaParsBreedFail)){
      warning("The combination of brood failure probability and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  if (varA > 0){
    betaParsA<- fbeta(mean=a, var=varA)
    if (anyNA(betaParsA)){
      warning("The combination of adult survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
    
    if (AFR > 1){
      betaParsS<- fbeta(mean=s, var=varA)
      if (anyNA(betaParsS)){
        warning("The combination of subadult survival and variance fall outside the domain of the Beta distribution.")
        return (NA)
      }
    }
  }
  
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  jEnv<- j
  breedFailEnv<- breedFail
  sEnv<- s
  aEnv<- a
  
  for (t in 1:tf){
    # Environmental stochasticity
    if (varJ > 0){
      jEnv<- rbeta(replicates, shape1=betaParsJ$shape1, shape2=betaParsJ$shape2)
    }
    if (varBreedFail > 0){
      breedFailEnv<- rbeta(replicates, shape1=betaParsBreedFail$shape1, shape2=betaParsBreedFail$shape2)
    }
    
    succeedingBroods<- rbinom(replicates, pop[, t, AFR] * broods, 1 - breedFailEnv)
    # Juvenile survivors
    pop[, t+1, 1]<- rbinom(replicates, b * succeedingBroods, jEnv)
    
    # Add adult survivors
    if (varA > 0){
      aEnv<- rbeta(replicates, shape1=betaParsA$shape1, shape2=betaParsA$shape2)
    }
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], aEnv)
    
    # Subadults (growth transitions)
    if (varA > 0 & AFR > 1)
      sEnv<- rbeta(replicates, shape1=betaParsS$shape1, shape2=betaParsS$shape2)
    
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, sEnv))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, sEnv))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=j)
mSurvBV.tseason<- function(broods, b, j, s, a, AFR, breedFail, seasonVar, N0, replicates, tf, maxN=100000){
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")

  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  for (t in 1:tf){
    for (k in 1:broods){
      succeedingBroods<- rbinom(replicates, pop[, t, AFR], jbrSeason[k])
      # Juvenile survivors
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, b * succeedingBroods, jindSeason[k])
    }
    
    # Add adult survivors
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], a)
    
    # Subadults (growth transitions)
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, s))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, s))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ))
mSurvBV.tvarseason<- function(broods, b, j, s, a, AFR, breedFail, seasonVar, varJ=0, varBreedFail=0, varA=0, N0, replicates, tf, maxN=100000){

  if (varJ > 0){
    betaParsJind<- fbeta(mean=j, var=varJ)
    if (anyNA(betaParsJind)){
      warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  if (varBreedFail > 0){
    betaParsBreedFail<- fbeta(mean=breedFail, var=varBreedFail)
    if (anyNA(betaParsBreedFail)){
      warning("The combination of brood failure probability and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  if (varA > 0){
    betaParsA<- fbeta(mean=a, var=varA)
    if (anyNA(betaParsA)){
      warning("The combination of adult survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
    
    if (AFR > 1){
      betaParsS<- fbeta(mean=s, var=varA)
      if (anyNA(betaParsS)){
        warning("The combination of subadult survival and variance fall outside the domain of the Beta distribution.")
        return (NA)
      }
    }
  }
  
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  
  jindEnv<- j
  breedFailEnv<- breedFail
  sEnv<- s
  aEnv<- a
  
  for (t in 1:tf){
    # Environmental stochasticity
    if (varJ > 0){
      jindEnv<- rbeta(replicates, shape1=betaParsJind$shape1, shape2=betaParsJind$shape2)
    }
    
    if (varBreedFail > 0){
      breedFailEnv<- rbeta(replicates, shape1=betaParsBreedFail$shape1, shape2=betaParsBreedFail$shape2)
    }
    
    for (k in 1:broods){
      succeedingBroods<- rbinom(replicates, pop[, t, AFR], (1 - breedFailEnv) * seasonVar[k])
      
      # Juvenile survivors
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, b * succeedingBroods, jindEnv * seasonVar[k])
    }
      
    # Add adult survivors
    if (varA > 0)
      aEnv<- rbeta(replicates, shape1=betaParsA$shape1, shape2=betaParsA$shape2)
    
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], aEnv)
    
    # Subadults (growth transitions)
    if (varA > 0 & AFR > 1)
      sEnv<- rbeta(replicates, shape1=betaParsS$shape1, shape2=betaParsS$shape2)
    
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, sEnv))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, sEnv))
    }
    
    pop[which(pop[, t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- cleanDiscretePopSim(pop, maxN=maxN)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=B(n=N_t, p=a) * fecundity, p=Beta(j, varJ)), p=sexRatio)
# fecundity = b x broods
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# TODO: add subadults
# TODO: varA
mFitSex.tvar<- function(fecundity, j, s, a, AFR, varJ=0, varA=0, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  ## TODO
  warning("Model not implemented?: mFitSex.tvar")

  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  
  betaPars<- fbeta(mean=j, var=varJ)
  if (anyNA(betaPars)){
    warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
    return (NA)
  }
  
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
# TODO: add subadults
# TODO: varA
mSurvBVSex.tvar<- function(broods, b, j, s, a, AFR, breedFail, varJ=0, varBreedFail=0, varA=0, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry"), N0, replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)
  
  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  if (varJ > 0){
    betaParsJ<- fbeta(mean=j, var=varJ)
    if (anyNA(betaParsJ)){
      warning("The combination of juvenile survival and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
  if (varBreedFail > 0){
    betaParsBreedFail<- fbeta(mean=breedFail, var=varBreedFail)
    if (anyNA(betaParsBreedFail)){
      warning("The combination of brood failure probability and variance fall outside the domain of the Beta distribution.")
      return (NA)
    }
  }
  
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
# TODO: add subadults
mSurvBVSex.tseason<- function(broods, b, j, s, a, AFR, breedFail, seasonVar, N0, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)

  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  popF<- popM<- matrix(NA_real_, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                    monogamy=min(popF[,t], popM[,t]),
                    polygyny=popF[,t],
                    polyandry=popM[,t])
    
    for (k in seq_along(broods)){
      succeedingBroods<- rbinom(replicates, nPairs, jbrSeason[k])
      # Juvenile survivors
      reclutes<- rbinom(1, succeedingBroods * b * 2, jindSeason[k])
      # Sex ratio
      popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
      popM[,t+1]<- reclutes - popF[,t+1]
    }
    
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
# N_t+1 = B(B(n=B(n=B(n=N_t, p=a) * broods, p=1-breedFail) * b, p=Beta(j, varJ)), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# TODO: add subadults
mSurvBVSex.tvarseason<- function(broods, b, j, s, a, AFR, breedFail, seasonVar, varJ=0, varBreedFail=0, varA=0, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), N0, replicates, tf, maxN=100000){
  matingSystem<- match.arg(matingSystem)
  
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  ## TODO
  warning("Model not implemented: mSurvBVSex.tvarseason")
  return (NA)
}

