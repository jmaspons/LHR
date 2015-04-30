## return object of class c("discretePopSim", "data.frame") with replicates on rows and time in columns

setOldClass("discretePopSim")

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
# N_t+1 = B(n=NB(n=N_t, p=1-a) * fecundity, p=j)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
mFit.t<- function(LH, N0, replicates, tf, maxN=100000){
  pop<- with(LH, expr={
    pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
    pop[,1]<- N0
    for (t in 1:tf){
      # Juvenile survivors
      pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, j)
      # Add adult survivors
      pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
      pop[which(pop[,t+1] > maxN),t+1]<- maxN
    }
    pop<- pop[order(pop[,tf+1]),]
    pop<- extinctNA(pop)
    pop<- as.data.frame(pop)
    class(pop)<- c("discretePopSim", "data.frame")
    pop
  })
  return(pop)
}

# Adult mortality + Nest mortality + offspring mortality
# N_t+1 = B(n=B(n=NB(n=N_t, p=1-a) * broods, p=1-nestFail) * b, p=j)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.t<- function(LH, nestFail, N0, replicates, tf, maxN=100000){
  LH$j<- LH$j / (1 - nestFail)
#   to-report annual2monthMortality [p]      
#   annualSurv = monthlySurv ^ 12    =>   monthlySurv = anualSurv ^ 1 / 12
#   monthSurv = (1 - annualMortality) ^ (1 / 12)
#   negative binomial: annualSurv = pnbinom(q=0, size=12, prob=monthlySurv) => 0 fails en 12 trails
#   report 1 - (1 - p) ^ (1 / 12)
  
  
  pop<- with(LH, expr={
    pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
    pop[,1]<- N0
    for (t in 1:tf){
      succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - nestFail)
      # Juvenile survivors
      pop[,t+1]<- rbinom(replicates, b * succeedingBroods, j)
      # Add adult survivors
      pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], a)
      pop[which(pop[,t+1] > maxN),t+1]<- maxN
    }
    pop<- pop[order(pop[,tf+1]),]
    pop<- extinctNA(pop)
    pop<- as.data.frame(pop)
    class(pop)<- c("discretePopSim", "data.frame")
    pop
  })
  return(pop)
}


# Adult mortality + offspring mortality + sex ratio 
# N_t+1 = B(n=B(n=NB(n=N_t, p=1-a) * fecundity, p=j), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = b x broods
## TODO:
mFitSex.t<- function(LH, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  pops<- with(LH, expr={
    popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
    popF[,1]<- popM[,1]<- N0
    for (t in 1:tf){
      nPairs<- switch(matingSystem, 
                      monogamy=min(popF[,t], popM[,t]),
                      polygyny=popF[,t],
                      polyandry=popM[,t])
      # Juvenile survivors
      reclutes<- rbinom(replicates, nPairs * fecundity, j)
      # Sex ratio
      popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
      popM[,t+1]<- reclutes - popF[,t+1]
      # Add adult survivors
      popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
      popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
      popF[which(popF[,t+1] > maxN),t+1]<- maxN
      popM[which(popM[,t+1] > maxN),t+1]<- maxN
    }
    popF<- extinctNA(popF)
    popM<- extinctNA(popM)
    popF<- as.data.frame(popF)
    popM<- as.data.frame(popM)
    class(popF)<- class(popM)<- c("discretePopSim", "data.frame")
    pops<- list(females=popF, males=popM)
  })
  return(pops)
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=NB(n=N_t, p=1-a) * broods, p=1-nestFail) * b, p=j), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.t<- function(LH, nestFail, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf, maxN=100000){
  LH$j<- LH$j / (1 - nestFail)
  pops<- with(LH, expr={
    popF<- popM<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
    popF[,1]<- popM[,1]<- N0
    for (t in 1:tf){
      nPairs<- switch(matingSystem, 
                       monogamy=min(popF[,t], popM[,t]),
                       polygyny=popF[,t],
                       polyandry=popM[,t])
      succeedingBroods<- rbinom(replicates, nPairs * broods, 1 - nestFail)
      # Juvenile survivors
      reclutes<- rbinom(replicates, b * succeedingBroods, j)
      # Sex ratio
      popF[,t+1]<- rbinom(replicates, reclutes, sexRatio)
      popM[,t+1]<- reclutes - popF[,t+1]
      # Add adult survivors
      popF[,t+1]<- popF[,t+1] + rbinom(replicates, popF[,t], a)
      popM[,t+1]<- popM[,t+1] + rbinom(replicates, popM[,t], a)
      popF[which(popF[,t+1] > maxN),t+1]<- maxN
      popM[which(popM[,t+1] > maxN),t+1]<- maxN
    }
    popF<- extinctNA(popF)
    popM<- extinctNA(popM)
    popF<- as.data.frame(popF)
    popM<- as.data.frame(popM)
    class(popF)<- class(popM)<- c("discretePopSim", "data.frame")
    pops<- list(females=popF, males=popM)
  })
  return(pops)
}


## SEASONAL ENVIRONMENT ---- TODO
