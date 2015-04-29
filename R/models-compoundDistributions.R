## Transition models N_t+1 = f(N_t) ----

# Adult mortality + offspring mortality
# N_t+1 = B(n=N_t * fecundity, p=j) + B(n=N_t, p=a)
mFit.trans<- function(LH, N){
  N_t1<- with(LH, expr={
    ## Recruitment probability conditioned to the number of nest success
    recruits<- distriBinom(N * fecundity, prob=j)
    
    ## Adult survival
    survivors<- distriBinom(size=N, prob=a)
  
    distriSum(recruits, survivors)
  })
  return (N_t1)
}

# Adult mortality + Nest mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t * broods, p=1-nestFail) * b, p=j) + B(n=N_t, p=a)
mSurvBV.trans<- function(LH, nestFail, N){
  LH$j<- nestFail / LH$j
  N_t1<- with(LH, expr={
    ## Probability of successful nesting without complete breeding fail
    nSuccessBroods<- with(LH, distriBinom(N * broods, prob=1-nestFail))
    nSuccessBroods$x<- nSuccessBroods$x * b
    ## Recruitment probability conditioned to the number of nest success
    recruits<- distriBinom(nSuccessBroods, prob=LH$j)
    
    ## Adult survival
    survivors<- distriBinom(size=N, prob=a)
    N_t1<- distriSum(recruits, survivors)
  })
  return (N_t1)
}

##TODO: what if there is differential mortalities between males and females?
## SEX RATIO: reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=N_t * fecundity, p=j), p=sexRatio) + B(n=N_t, p=a)
mFitSex.trans<- function(LH, sexRatio=0.5, Nf, Nm, matingSystem=c("monogamy", "polygyny", "polyandry")[1]){
  N_t1<- with(LH, expr={
    ## Recruitment probability conditioned to the number of nest success
    N<- switch(matingSystem, 
                    monogamy=min(Nf, Nm),
                    polygyny=Nf,
                    polyandry=Nm)
    recruits<- distriBinom(N * fecundity, prob=j)
    recruitsF<- distriBinom(recruits, prob=sexRatio)
    recruitsM<- recruits - recruitsF
    
    ## Adult survival
    survivorsF<- distriBinom(size=Nf, prob=a)
    survivorsM<- distriBinom(size=NM, prob=a)
    Nf_t1<- distriSum(recruitsF, survivorsF)
    Nm_t1<- distriSum(recruitsF, survivorsM)
  })
  return (N_t1)
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=N_t * broods, p=1-nestFail) * b, p=j), p=sexRatio) + B(n=N_t, p=a)
mSurvBVSex.trans<- function(LH, nestFail, sexRatio=0.5, Nf, Nm){
  LH$j<- nestFail / LH$j
  N_t1<- with(LH, expr={  
    
  })
}


## Cohort models ----
# Expected life long fitness for a cohort.

# Nest mortality + offspring mortality
# cohort fitness = B(B(n=broods * lifespan, p=1-nestFail) * b, p=jusSurv)
mBV.cohort<- function(LH, lifespan, nestFail){
  LH$j<- nestFail / LH$j
  N_t1<- with(LH, expr={
    ## Probability of successful nesting without complete breeding fail
    nSuccessBroods<- distriBinom(broods * lifespan, prob=1-nestFail)
    # Multiply the number of successful broods by b size
    nSuccessBroods$x<- nSuccessBroods$x * b
    ## Fitness probability conditioned to the number of nest success
    fitness<- distriBinom(nSuccessBroods, prob=j)
  })
  return (fitness)
}

# Adult mortality + offspring mortality
# cohort fitness = B(n=NB(n=N_0, p=1-a) * fecundity, p=j)
mFit.cohort<- function(LH, N0=1, ...){ # ... = maxX or p.omited 
  N_t1<- with(LH, expr={
    lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
    fitness<- distriBinom(size=fecundity, prob=j)
  })
  return (fitness)
}

# Adult mortality + Nest mortality + offspring mortality
# cohort fitness = B(n=B(n=NB(n=N_0, p=1-a) * broods, p=1-nestFail) * b, p=j)
mSurvBV.cohort<- function(LH, nestFail, N0=1, ...){ # ... = maxX or p.omited
  LH$j<- nestFail / LH$j
  N_t1<- with(LH, expr={
    lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
    lifespan$x<- lifespan$x * broods # Number of broods for the whole life
    ## Probability of successful nesting without complete breeding fail
    nSuccessBroods<- distriBinom(lifespan, prob=1-nestFail)
    nSuccessBroods$x<- nSuccessBroods$x * b
    ## Fitness probability conditioned to the number of nest success
    fitness<- distriBinom(nSuccessBroods, prob=j)
  })
  return (fitness)
}

##TODO: how to deal with sex ratio??? Number of males and females are not independent, they sum fitness and only change the proportion.
## Reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=NB(n=N_0, p=1-a) * fecundity, p=j), p=sexRatio)
mFitSex.cohort<- function(LH, sexRatio=0.5, N0=1, ...){ # ... = maxX or p.omited 
  N_t1<- with(LH, expr={
    lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
    fitness<- distriBinom(size=fecundity, prob=j)
    females<- distriBinom(size=fitness, prob=sexRatio)
    males<- distriBinom(size=fitness, prob=1-sexRatio)
  })
  return (list(females=females, males=males))
}

# Adult mortality + Nest mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=B(n=NB(n=N_0, p=1-a) * broods, p=1-nestFail) * b, p=j), p=sexRatio)
mSurvBVSex.cohort<- function(LH, nestFail, sexRatio=0.5, N0=1, ...){ # ... = maxX or p.omited 
  LH$j<- nestFail / LH$j
  N_t1<- with(LH, expr={
    lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
    lifespan$x<- lifespan$x * broods # Number of broods for the whole life
    ## Probability of successful nesting without complete breeding fail
    nSuccessBroods<- distriBinom(lifespan, prob=1-nestFail)
    nSuccessBroods$x<- nSuccessBroods$x * b
    ## Fitness probability conditioned to the number of nest success
    fitness<- distriBinom(nSuccessBroods, prob=j)
    
    females<- distriBinom(size=fitness, prob=sexRatio)
    males<- distriBinom(size=fitness, prob=1-sexRatio)
  })
  return (list(females=females, males=males))
}

