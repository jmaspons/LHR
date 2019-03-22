## Multiple time steps models decrease accuracy after 50 timesteps aprox. Use logP(p) reduce omitted probability
## TODO: fix log probabilities when = 1 or 0??

#' @include numericDistri.R
NULL

## S4 model dispatcher ----

#' Numeric distribution models
#' 
#' Calculate the distribution of the population size for a transition of \code{tf} time steps.
#'
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
#' @param tf 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("tDistri", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, seasonVar=0, sexRatio=NA,
                                matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), N0=10, tf=1) standardGeneric("tDistri"))
setMethod("tDistri",
          signature(broods="ANY", b="numeric", j="numeric", a="numeric", breedFail="ANY", varJ="ANY", varBreedFail="ANY",
                    seasonVar="ANY", sexRatio="ANY", matingSystem="ANY", N0="ANY", tf="ANY"),
          function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, seasonVar=0, sexRatio=.5,
                   matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), N0=10, tf=1){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }

            out<- NA
            
            tf1<- broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE

            if (tf == 1) tf1<- TRUE                               # only one time step
            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0) varEnv<- TRUE       # environmental variation
            if (length(j) > 1 | length(a) > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes
            
            
            if (!broodMortality & !varEnv & !seasonalEnv & !sex){
              out<- mFit.distri(fecundity=b * broods, j=j, a=a, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mFit.distri(fecundity=b * broods, j=j, a=a, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
              
            }
            if (!broodMortality & !varEnv & !seasonalEnv &  sex){
              out<- mFitSex.distri(fecundity=b * broods, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mFitSex.distri(fecundity=b * broods, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, Nf=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (!broodMortality & !varEnv & seasonalEnv  & !sex){
              out<- mFit_season.distri(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mFit_season.distri(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (!broodMortality & !varEnv & seasonalEnv  &  sex){
              ## TODO 
              warning("## TODO: out<- mFitSex_season.distri")
            }
            if (!broodMortality & varEnv  & !seasonalEnv & !sex){
              out<- mFit_var.distri(fecundity=b * broods, j=j, a=a, varJ=varJ, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mFit_var.distri(fecundity=b * broods, j=j, a=a, varJ=varJ, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (!broodMortality & varEnv  & !seasonalEnv &  sex){
              ## TODO 
              warning("## TODO: out<- mFitSex_var.distri")
              # out<- mFitSex_var.distri(fecundity=b * broods, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
            }
            if (!broodMortality & varEnv  & seasonalEnv  & !sex){
              out<- mFit_varseason.distri(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, varJ=varJ, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mFit_varseason.distri(broods=broods, b=b, j=j, a=a, seasonVar=seasonVar, varJ=varJ, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (!broodMortality & varEnv  & seasonalEnv  &  sex){
              ## TODO 
              warning("## TODO: out<- mFitSex_varseason.distri")
            }
            
            if (broodMortality  & !varEnv & !seasonalEnv & !sex){
              out<- mSurvBV.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mSurvBV.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (broodMortality  & !varEnv & !seasonalEnv &  sex){
              out<- mSurvBVSex.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mSurvBVSex.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, Nf=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (broodMortality  & !varEnv & seasonalEnv  & !sex){
              out<- mSurvBV_season.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mSurvBV_season.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (broodMortality  & !varEnv & seasonalEnv  &  sex){
              ## TODO 
              warning("## TODO: out<- mSurvBVSex_season.distri")
              # out<- mSurvBVSex_season.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
            }
            if (broodMortality  & varEnv  & !seasonalEnv & !sex){
              out<- mSurvBV_var.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mSurvBV_var.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (broodMortality  & varEnv  & !seasonalEnv &  sex){
              ## TODO 
              warning("## TODO: out<- mSurvBVSex_var.distri")
              # out<- mSurvBVSex_var.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
            }
            if (broodMortality  & varEnv  & seasonalEnv  & !sex){
              out<- mSurvBV_varseason.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, N0=N0)
              at<- attributes(out)
              
              if (tf > 1){
                for (i in 2:tf){
                  out<- mSurvBV_varseason.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, seasonVar=seasonVar, varJ=varJ, varBreedFail=varBreedFail, N0=out)
                  out<- out[!out$p %in% c(-Inf, 0),]
                }
              }
            }
            if (broodMortality  & varEnv  & seasonalEnv  &  sex){
              ## TODO 
              warning("## TODO: out<- mSurvBVSex_varseason.distri")
              # out<- mSurvBVSex_varseason.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, Nf=N0)
            }
            
            attributes(out)$p.omitted<- 1 - sum(logP(out, logP=FALSE)$p)
            attributes(out)$parameters<- at$parameters
            attributes(out)$tf<- tf
            attributes(out)$N0<- N0
            
            class(out)<- c("numericDistriPopSim", class(out))
            
            return (out)
          }
)


## Transition models N_t+1 = f(N_t) ----
# Adult mortality + offspring mortality
# N_t+1 = B(n=N_t * fecundity, p=j) + B(n=N_t, p=a)
mFit.distri<- function(fecundity, j, a, N0, logP=FALSE){
  ## Recruitment probability conditioned to the number of successful broods
  recruits<- distriBinom(size=N0 * fecundity, prob=j, logP=logP)
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors

  return (N_t1)
}

# N_t+1 = B(n=N_t * fecundity, p=Beta(j, varJ)) + B(n=N_t, p=a)
mFit_var.distri<- function(fecundity, j, a, varJ, N0, logP=FALSE){
  betaPars<- fbeta(mean=j, var=varJ)
  # if (all(is.na(betaPars)))
  recruits<- distriBetaBinom(N0 * fecundity, shape1=betaPars$shape1, shape2=betaPars$shape2, logP=logP)
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors
  
  return (N_t1)
}


mFit_season.distri<- function(broods, b, j, a, seasonVar, N0, logP=FALSE){
  jindSeason<- j * seasonVar
  
  recruits<- 0
  for (k in 1:broods){
    # Juvenile survivors
    recruits<- recruits + distriBinom(size=b , prob=jindSeason[k], logP=logP)
  }
  
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors

  return (N_t1)
}

mFit_varseason.distri<- function(broods, b, j, a, seasonVar, varJ, N0){
  betaJ<- fbeta(mean=j, var=varJ)

  recruits<- 0
  for (k in 1:broods){
    # Juvenile survivors
    recruits<- recruits + distriBetaBinom(size=b , shape1=betaJ$shape1, shape2=betaJ$shape2, logP=logP)
  }

  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors

  return (N_t1)
}


## Brood Mortality ----
# Adult mortality + Brood mortality + offspring mortality
# N_t+1 = B(n=B(n=N_t * broods, p=1-breedFail) * b, p=j) + B(n=N_t, p=a)
mSurvBV.distri<- function(broods, b, j, a, breedFail, N0, logP=FALSE){
  ## Probability of successful breeding attempt without complete brood fail
  nSuccessBroods<- distriBinom(N0 * broods, prob=1-breedFail, logP=logP)
  nSuccessBroods<- nSuccessBroods * b
  
  ## Recruitment probability conditioned to the number of successful broods
  recruits<- distriBinom(nSuccessBroods, prob=j, logP=logP)
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors
  
  return (N_t1)
}

# N_t+1 = B(n=B(n=N_t * broods, p=Beta(1-breedFail, varBreedFail)) * b, p=Beta(j, varJ)) + B(n=N_t, p=a)
mSurvBV_var.distri<- function(broods, b, j, a, breedFail, varJ, varBreedFail, N0, logP=FALSE){
  ## Probability of successful breeding attempt without complete brood fail
  if (varBreedFail > 0){
    betaPars<- fbeta(mean=1-breedFail, var=varBreedFail)
    nSuccessBroods<- distriBetaBinom(N0 * broods, shape1=betaPars$shape1, shape2=betaPars$shape2, logP=logP)
  }else{
    nSuccessBroods<- distriBinom(N0 * broods, prob=1-breedFail, logP=logP)
  }
  
  nSuccessBroods<- nSuccessBroods * b
  
  ## Recruitment probability conditioned to the number of successful broods
  if (varJ > 0){
    betaPars<- fbeta(mean=j, var=varJ)
    recruits<- distriBetaBinom(nSuccessBroods, shape1=betaPars$shape1, shape2=betaPars$shape2, logP=logP)
  }else{
    recruits<- distriBinom(nSuccessBroods, prob=j, logP=logP)
  }
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors
  
  return (N_t1)
}

# N_t+1 = B(n=B(n=N_t * broods, p=Beta(1-breedFail, varBreedFail)) * b, p=Beta(j, varJ)) + B(n=N_t, p=a)
mSurvBV_season.distri<- function(broods, b, j, a, breedFail, seasonVar, N0, logP=FALSE){
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  recruits<- 0
  for (k in 1:broods){
    ## Probability of successful breeding attempt without complete brood fail
    nSuccessBroods<- distriBinom(size=N0, prob=jbrSeason[k], logP=logP)
    nSuccessBroods<- nSuccessBroods * b

    ## Recruitment probability conditioned to the number of successful broods
    recruits<- distriBinom(nSuccessBroods, prob=jindSeason[k], logP=logP)
  }
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors
  
  return (N_t1)
}

# N_t+1 = B(n=B(n=N_t * broods, p=Beta(1-breedFail, varBreedFail)) * b, p=Beta(j, varJ)) + B(n=N_t, p=a)
mSurvBV_varseason.distri<- function(broods, b, j, a, breedFail, seasonVar, varJ, varBreedFail, N0, logP=FALSE){
  jindSeason<- j * seasonVar
  jbrSeason<- (1 - breedFail) * seasonVar
  
  betaJind<- fbeta(mean=jindSeason, var=varJ)
  betaJbr<- fbeta(mean=jbrSeason, var=varBreedFail)
  
  recruits<- 0
  for (k in 1:broods){
    ## Probability of successful breeding attempt without complete brood fail
    nSuccessBroods<- distriBetaBinom(size=N0, shape1=betaJbr$shape1[k], shape2=betaJbr$shape2[k], logP=logP)
    nSuccessBroods<- nSuccessBroods * b

    ## Recruitment probability conditioned to the number of successful broods
    recruits<- distriBetaBinom(nSuccessBroods, shape1=betaJind$shape1[k], shape2=betaJind$shape2[k], logP=logP)
  }
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- recruits + survivors
  
  return (N_t1)
}


## 2 sexes ----
##TODO: what if there is differential mortalities between males and females?
## SEX RATIO: reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=N_t * fecundity, p=j), p=sexRatio) + B(n=N_t, p=a)
mFitSex.distri<- function(fecundity, j, a, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), Nf, Nm=Nf, logP=FALSE){
  matingSystem<- match.arg(matingSystem)
  
  ## Recruitment probability conditioned to the number of successful broods
  N<- switch(matingSystem, 
                  monogamy=min(Nf, Nm),
                  polygyny=Nf,
                  polyandry=Nm)
  recruits<- distriBinom(N * fecundity, prob=j, logP=logP)
  recruitsF<- distriBinom(recruits, prob=sexRatio, logP=logP)
  recruitsM<- recruits - recruitsF
  
  ## Adult survival
  survivorsF<- distriBinom(size=Nf, prob=a, logP=logP)
  survivorsM<- distriBinom(size=Nm, prob=a, logP=logP)
  Nf_t1<- distriSum(recruitsF, survivorsF)
  Nm_t1<- distriSum(recruitsF, survivorsM)
  N_t1<- list(females=Nf_t1, males=Nm_t1)
    
  return (N_t1)
}

# Adult mortality + Brood mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=B(n=N_t * broods, p=1-breedFail) * b, p=j), p=sexRatio) + B(n=N_t, p=a)
mSurvBVSex.distri<- function(broods, b, j, a, breedFail, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry"), Nf, Nm, logP=FALSE){
  matingSystem<- match.arg(matingSystem)

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

#' @export
lambda.numericDistri<- function(x, N0, tf=1, ...){
  x$x<- (x$x / N0)^(1/tf)
  
  if (N0 == 0){
    x$x<- 0
    x$p<- 1
  }
  
  return(x)
}

#' @export
r.numericDistri<- function(x, N0, tf=1, ...){
  x$x<- (x$x - N0) / N0 / tf
  
  if (N0 == 0){
    x$x<- -Inf
    x$p<- 1
  }
  
  return(x)
}
