# Nest mortality + offspring mortality
# fitness = B(B(n=broods * lifespan, p=1-nestFail) * clutch, p=jusSurv)
mBV.dist<- function(broods, clutch, lifespan, nestFail, juvSurv){
  ## Probability of successful nesting without complete breeding fail
  nSuccessBroods<- distriBinom(broods * lifespan, prob=1-nestFail)
  # Multiply the number of successful broods by clutch size
  nSuccessBroods$x<- nSuccessBroods$x * clutch
  ## Fitness probability conditioned to the number of nest success
  fitness<- distriBinom(nSuccessBroods, prob=juvSurv)
  
  return (fitness)
}

# Adult mortality + offspring mortality
# fitness = B(n=NB(n=N_0, p=1-adultSurv) * fecundity, p=juvSurv)
mFit.dist<- function(fecundity, juvSurv, adultSurv, N0=1, ...){ # ... = maxX or p.omited 
  lifespan<- distriNegBinom(size=N0, prob=1-adultSurv, ...)
  fitness<- distriBinom(size=fecundity, prob=juvSurv)
  
  return (fitness)
}


# Adult mortality + Nest mortality + offspring mortality
# fitness = B(n=B(n=NB(n=N_0, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv)
mSurvBV.dist<- function(broods, clutch, nestFail, juvSurv, adultSurv, N0=1, ...){ # ... = maxX or p.omited 
  lifespan<- distriNegBinom(size=N0, prob=1-adultSurv, ...)
  lifespan$x<- lifespan$x * broods # Number of broods for the whole life
  ## Probability of successful nesting without complete breeding fail
  nSuccessBroods<- distriBinom(lifespan, prob=1-nestFail)
  nSuccessBroods$x<- nSuccessBroods$x * clutch
  ## Fitness probability conditioned to the number of nest success
  fitness<- distriBinom(nSuccessBroods, prob=juvSurv)
  
  return (fitness)
}

##TODO: how to deal with sex ratio??? Number of males and females are not independent, they sum fitness and only change the proportion.
## Reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=NB(n=N_0, p=1-adultSurv) * fecundity, p=juvSurv), p=sexRatio)
mFitSex.dist<- function(fecundity, juvSurv, adultSurv, sexRatio=0.5, N0=1, ...){ # ... = maxX or p.omited 
  lifespan<- distriNegBinom(size=N0, prob=1-adultSurv, ...)
  fitness<- distriBinom(size=fecundity, prob=juvSurv)
  females<- distriBinom(size=fitness, prob=sexRatio)
  males<- distriBinom(size=fitness, prob=1-sexRatio)
  return (list(females=females, males=males))
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=B(n=NB(n=N_0, p=1-adultSurv) * broods, p=1-nestFail) * clutch, p=juvSurv), p=sexRatio)
mSurvBVSex.dist<- function(broods, clutch, nestFail, juvSurv, adultSurv, sexRatio=0.5, N0=1, ...){ # ... = maxX or p.omited 
  lifespan<- distriNegBinom(size=N0, prob=1-adultSurv, ...)
  lifespan$x<- lifespan$x * broods # Number of broods for the whole life
  ## Probability of successful nesting without complete breeding fail
  nSuccessBroods<- distriBinom(lifespan, prob=1-nestFail)
  nSuccessBroods$x<- nSuccessBroods$x * clutch
  ## Fitness probability conditioned to the number of nest success
  fitness<- distriBinom(nSuccessBroods, prob=juvSurv)
  
  females<- distriBinom(size=fitness, prob=sexRatio)
  males<- distriBinom(size=fitness, prob=1-sexRatio)
  return (list(females=females, males=males))
}
