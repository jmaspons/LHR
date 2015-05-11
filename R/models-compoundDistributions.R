## S4 model dispatcher ----
setGeneric("t1distri", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, sexRatio=.5,
                                matingSystem=c("monogamy", "polygyny", "polyandry")[1], N) standardGeneric("t1distri"))

setMethod("t1distri",  # function dispatcher # TODO check what happen with the full model and the dispatcher (same signature)
          signature(broods="numeric", b="numeric", j="numeric", a="numeric", breedFail="numeric", varJ="numeric", varBreedFail="numeric",
                    sexRatio="ANY", matingSystem="ANY", N="numeric"),
          function(broods=1, b, j, a, breedFail, varJ=0, varBreedFail=0, sexRatio=.5, matingSystem="monogamy", N){
            cm<- paste0("t1distri(j=j, a=a, b=b, ")
            if (any(breedFail != 0)){ # any for seasonal where j and breedFail can be a vector of broods length
              cm<- paste0(cm, "broods=broods, breedFail=breedFail, ")
            }
            if (!is.na(sexRatio)){
              cm<- paste0(cm, "sexRatio=sexRatio, matingSystem=matingSystem, ")
            }
            if (varJ != 0 | varBreedFail !=0){
              # Select one of the following options and think about correlation (same mean for juvenile survival and breeding fail?)
              cm<- paste0(cm, "varJ=varJ, varBreedFail=0, ")
              #               cm<- paste0(cm, "varJ=0, varBreedFail=varBreedFail, ")
              #               cm<- paste0(cm, "varJ=varJ, varBreedFail=varBreedFail,")
            }
            
            cm<- paste0(cm, "N=N)")
            print(cm)
            eval(expr=parse(text=cm))
          }
)

## Stable environment
setMethod("t1distri", 
          signature(broods="missing", b="numeric", j="numeric", a="numeric", breedFail="missing", varJ="missing", varBreedFail="missing",
                    sexRatio="missing", matingSystem="missing", N="numeric"),
          function(b, j, a, N){ #mFit
            mFit.trans(fecundity=b, j=j, a=a, N=N)
          }
)
setMethod("t1distri", 
          signature(broods="numeric", b="numeric", j="numeric", a="numeric", breedFail="numeric", varJ="missing", varBreedFail="missing",
                    sexRatio="missing", matingSystem="missing", N="numeric"),
          function(broods, b, j, a, breedFail, N){ #mSurvBV
            if (length(j) > 1 | length(a) > 1){
              mSurvBV.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N=N)
            }else{
              mSurvBV.trans(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N=N)
            }
          }
)
setMethod("t1distri", 
          signature(broods="missing", b="numeric", j="numeric", a="numeric", breedFail="missing", varJ="missing", varBreedFail="missing",
                    sexRatio="numeric", matingSystem="character", N="numeric"),
          function(b, j, a, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N){ #mFitSex
            mFitSex.trans(fecundity=b, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N=N)
          }
)
setMethod("t1distri", 
          signature(broods="numeric", b="numeric", j="numeric", a="numeric", breedFail="numeric", varJ="missing", varBreedFail="missing",
                    sexRatio="numeric", matingSystem="character", N="numeric"),
          function(broods, b, j, a, breedFail, sexRatio=.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], N){ #mSurvBVSex
            if (length(j) > 1 | length(a) > 1){
              mSurvBVSex.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N=N)
            }else{
              mSurvBVSex.trans(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N=N)
            }
          }
)

# Stochastic environment
#TODO


## Transition models N_t+1 = f(N_t) ----
# Adult mortality + offspring mortality
# N_t+1 = B(n=N_t * fecundity, p=j) + B(n=N_t, p=a)
mFit.trans<- function(fecundity, j, a, N){
  ## Recruitment probability conditioned to the number of successful broods
  recruits<- distriBinom(N * fecundity, prob=j)
  
  ## Adult survival
  survivors<- distriBinom(size=N, prob=a)
  N_t1<- distriSum(recruits, survivors)

  return (N_t1)
}

# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t * broods, p=1-breedFail) * b, p=j) + B(n=N_t, p=a)
mSurvBV.trans<- function(broods, b, j, a, breedFail, N){
  ## Probability of successful breeding attempt without complete brood fail
  nSuccessBroods<- distriBinom(N * broods, prob=1-breedFail)
  nSuccessBroods$x<- nSuccessBroods$x * b
  ## Recruitment probability conditioned to the number of successful broods
  recruits<- distriBinom(nSuccessBroods, prob=j)
  
  ## Adult survival
  survivors<- distriBinom(size=N, prob=a)
  N_t1<- distriSum(recruits, survivors)

  return (N_t1)
}

##TODO: what if there is differential mortalities between males and females?
## SEX RATIO: reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=N_t * fecundity, p=j), p=sexRatio) + B(n=N_t, p=a)
mFitSex.trans<- function(fecundity, j, a, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], Nf, Nm=Nf){
  ## Recruitment probability conditioned to the number of successful broods
  N<- switch(matingSystem, 
                  monogamy=min(Nf, Nm),
                  polygyny=Nf,
                  polyandry=Nm)
  recruits<- distriBinom(N * fecundity, prob=j)
  recruitsF<- distriBinom(recruits, prob=sexRatio)
  recruitsM<- recruits - recruitsF
  
  ## Adult survival
  survivorsF<- distriBinom(size=Nf, prob=a)
  survivorsM<- distriBinom(size=Nm, prob=a)
  Nf_t1<- distriSum(recruitsF, survivorsF)
  Nm_t1<- distriSum(recruitsF, survivorsM)
  N_t1<- list(females=Nf_t1, males=Nm_t1)
    
  return (N_t1)
}

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=N_t * broods, p=1-breedFail) * b, p=j), p=sexRatio) + B(n=N_t, p=a)
mSurvBVSex.trans<- function(broods, b, j, a, breedFail, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], Nf, Nm){

}


## Cohort models ----
# Expected life long fitness for a cohort.

# Brood mortality + offspring mortality
# cohort fitness = B(B(n=broods * lifespan, p=1-breedFail) * b, p=jusSurv)
mBV.cohort<- function(j, lifespan, breedFail){
  ## Probability of successful breeding attempt without complete brood fail
  nSuccessBroods<- distriBinom(broods * lifespan, prob=1-breedFail)
  # Multiply the number of successful broods by b size
  nSuccessBroods$x<- nSuccessBroods$x * b
  ## Fitness probability conditioned to the number of successful broods
  fitness<- distriBinom(nSuccessBroods, prob=j)

  return (fitness)
}

# Adult mortality + offspring mortality
# cohort fitness = B(n=NB(n=N_0, p=1-a) * fecundity, p=j)
mFit.cohort<- function(fecundity, j, a, N0=1, ...){ # ... = maxX or p.omited 
  lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
  fitness<- distriBinom(size=fecundity, prob=j)

  return (fitness)
}

# Adult mortality + Brood mortality + offspring mortality
# cohort fitness = B(n=B(n=NB(n=N_0, p=1-a) * broods, p=1-breedFail) * b, p=j)
mSurvBV.cohort<- function(broods, b, j, a, breedFail, N0=1, ...){ # ... = maxX or p.omited
  lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
  lifespan$x<- lifespan$x * broods # Number of broods for the whole life
  ## Probability of successful breeding attempt without complete brood fail
  nSuccessBroods<- distriBinom(lifespan, prob=1-breedFail)
  nSuccessBroods$x<- nSuccessBroods$x * b
  ## Fitness probability conditioned to the number of successful broods
  fitness<- distriBinom(nSuccessBroods, prob=j)

  return (fitness)
}

##TODO: how to deal with sex ratio??? Number of males and females are not independent, they sum fitness and only change the proportion.
## Reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=NB(n=N_0, p=1-a) * fecundity, p=j), p=sexRatio)
# mFitSex.cohort<- function(fecundity, j, a, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], Nf, Nm=Nf, ...){ # ... = maxX or p.omited 
#     lifespan<- distriNegBinom(size=N0, prob=1-a, ...)
#     fitness<- distriBinom(size=fecundity, prob=j)
#     females<- distriBinom(size=fitness, prob=sexRatio)
#     males<- distriBinom(size=fitness, prob=1-sexRatio)
# 
#   return (list(females=females, males=males))
# }

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# fitness = B(n=B(n=B(n=NB(n=N_0, p=1-a) * broods, p=1-breedFail) * b, p=j), p=sexRatio)


## Distribution stats ----
lambda.numericDistri<- function(distri, N0){
  distriLambda<- distri
  distriLambda$x<- distri$x / N0
  return(distriLambda)
}

r.numericDistri<- function(distri, N0){
  distriR<- distri
  distriR$x<- (distri$x - N0) / N0
  return(distriR)
}
