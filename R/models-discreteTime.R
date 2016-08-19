
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
## TODO
mSurvBVSex.tvarseason<- function(broods, b, j, a, breedFail, varJ=0, varBreedFail=0, N0, replicates, tf, maxN=100000){
  
}
