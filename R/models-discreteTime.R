## return object of class c("discretePopSim", "data.frame") with replicates on rows and time in columns
# When the population gets extinct it fills results with NAs.
extinctNA<- function(pop){
  extT<- apply(pop, 1, function(x) match(0, x))
  extT<- extT[extT < ncol(pop)]
  extPop<- which(!is.na(extT))
  for (i in seq_along(extPop)){
    pop[i,(extT[i]+1):ncol(pop)]<- NA
  }
  return(pop)
}

##  CONSTANT ENVIRONMENT
# Adult mortality + offspring mortality
# N_t+1 = B(n=NB(n=N_t, p=1-adultSurv) * fecundity, p=juvSurv)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = clutch x broods
mFit.t<- function(fecundity, juvSurv, adultSurv, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,tf+1]),]
  pop<- extinctNA(pop)
  pop<- as.data.frame(pop)
  class(pop)<- c("discretePopSim", "data.frame")
  
  return(pop)
}

# Adult mortality + Nest mortality + offspring mortality
# N_t+1 = B(n=B(n=NB(n=N_t, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.t<- function(broods, clutch, nestFail, juvSurv, adultSurv, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - nestFail)
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, clutch * succeedingBroods, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,tf+1]),]
  pop<- extinctNA(pop)
  pop<- as.data.frame(pop)
  class(pop)<- c("discretePopSim", "data.frame")
  
  return(pop)
}


# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=NB(n=N_t, p=1-adultSurv) * fecundity, p=juvSurv), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = clutch x broods
## TODO:
mFitSex.t<- function(fecundity, juvSurv, adultSurv, sexRatio=.5, N0, replicates, tf, maxN=100000){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
    pop[which(pop[,t+1] > maxN),t+1]<- maxN
  }
  pop<- pop[order(pop[,tf+1]),]
  pop<- extinctNA(pop)
  pop<- as.data.frame(pop)
  class(pop)<- c("discretePopSim", "data.frame")
  return(pop)
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=NB(n=N_t, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.t<- function(broods, clutch, nestFail, juvSurv, adultSurv, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  popF[,1]<- popM[,1]<- N0
  for (t in 1:tf){
    nPairs<- switch(matingSystem, 
                     monogamy=min(popF[,t], popM[,t]),
                     polygyny=popF[,t],
                     polyandry=popM[,t])
    succeedingBroods<- rbinom(replicates, nPairs * broods, 1 - nestFail)
    # Juvenile survivors
    reclutes<- rbinom(replicates, clutch * succeedingBroods, juvSurv)
    # Sex ratio
    popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
    popM[,t+1]<- reclutes - popF[,t+1]
    # Add adult survivors
    popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], adultSurv)
    popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], adultSurv)
    popF[which(popF[,t+1] > maxN),t+1]<- maxN
    popM[which(popM[,t+1] > maxN),t+1]<- maxN
  }
  popF<- extinctNA(popF)
  popM<- extinctNA(popM)
  popF<- as.data.frame(popF)
  popM<- as.data.frame(popM)
  class(popF)<- class(popM)<- c("discretePopSim", "data.frame")
  return(list(females=popF, males=popM))
}


## SEASONAL ENVIRONMENT ---- TODO
