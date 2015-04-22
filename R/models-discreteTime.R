## return object of class c("discretePopSim", "matrix") with replicates on rows and time in columns

# TODO: equivalent model with fixed lifespan requires to know the age of individuals
# # Nest mortality + offspring mortality
# # fitness = B(B(n=broods * lifespan, p=1-nestFail) * clutch, p=jusSurv)
# # Juveniles reach adult stage in one time step (age at first reproduction = 1)
# mBV.t<- function(broods, clutch, nestFail, adultSurv, juvSurv, N0, replicates, tf){
#   pop<- matrix(NA, replicates, tf, dimnames=list(replicate=NULL, t=NULL))
#   pop[,1]<- N0
#   for (t in 1:tf){
#     # breeding success
#     succeedingbroods<- rbinom(replicates, pop[,t] * broods, 1 - nestFail)
#     # Juvenile survival
#     pop[,t+1]<- rbinom(replicates, clutch * succeedingbroods, juvSurv)
#     # Add adult survivors
#     pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
#   }
#   return(pop[order(pop[,tf]),])
# }

# Adult mortality + offspring mortality
# fitness = B(n=NB(n=N_0, p=1-adultSurv) * fecundity, p=juvSurv)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = clutch x broods
mFit.t<- function(fecundity, juvSurv, adultSurv, N0, replicates, tf){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
  }
  pop[order(pop[,tf+1]),]
  class(pop)<- c("discretePopSim", "matrix")
  return(pop)
}

# Adult mortality + Nest mortality + offspring mortality
# fitness = B(n=B(n=NB(n=N_0, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBV.t<- function(broods, clutch, nestFail, juvSurv, adultSurv, N0, replicates, tf){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    succeedingBroods<- rbinom(replicates, pop[,t] * broods, 1 - nestFail)
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, clutch * succeedingBroods, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
  }
  pop[order(pop[,tf+1]),]
  class(pop)<- c("discretePopSim", "matrix")
  return(pop)
}


# Adult mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=NB(n=N_0, p=1-adultSurv) * fecundity, p=juvSurv), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
# fecundity = clutch x broods
mFitSex.t<- function(fecundity, juvSurv, adultSurv, sexRatio=.5, N0, replicates, tf){
  pop<- matrix(NA, replicates, tf+1, dimnames=list(replicate=NULL, t=0:tf))
  pop[,1]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[,t+1]<- rbinom(replicates, pop[,t] * fecundity, juvSurv)
    # Add adult survivors
    pop[,t+1]<- pop[,t+1] + rbinom(replicates, pop[,t], adultSurv)
  }
  pop[order(pop[,tf+1]),]
  class(pop)<- c("discretePopSim", "matrix")
  return(pop)
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=B(n=NB(n=N_0, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv), p=sexRatio)
# Juveniles reach adult stage in one time step (age at first reproduction = 1)
mSurvBVSex.t<- function(broods, clutch, nestFail, juvSurv, adultSurv, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N0, replicates, tf){
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
  }
  class(popF)<- class(popM)<- c("discretePopSim", "matrix")
  return(list(females=popF, males=popM))
}
