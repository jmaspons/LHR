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
#' @param sexRatio 
#' @param matingSystem 
#' @param N0 
#' @param tf 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("tDistri", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, sexRatio=NA,
                                matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), N0=10, tf=1) standardGeneric("tDistri"))
setMethod("tDistri",
          signature(broods="ANY", b="numeric", j="numeric", a="numeric", breedFail="ANY", varJ="ANY", varBreedFail="ANY",
                    sexRatio="ANY", matingSystem="ANY", N0="ANY", tf="ANY"),
          function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, sexRatio=.5,
                   matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), N0=10, tf=1){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }

            out<- NA
            
            broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE

            if (tf == 1) tf1<- TRUE                               # only one time step
            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0) varEnv<- TRUE       # environmental variation
            if (length(j) > 1 | length(a) > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes

            if (tf1){
              if (!broodMortality & !varEnv & !seasonalEnv & !sex){
                out<- mFit.distri(fecundity=b, j=j, a=a, N0=N0)
              }
              if (!broodMortality & !varEnv & !seasonalEnv &  sex){
                out<- mFitSex.distri(fecundity=b, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              if (!broodMortality & !varEnv & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mFit_season.distri")
              }
              if (!broodMortality & !varEnv & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_season.distri")
              }
              if (!broodMortality & varEnv  & !seasonalEnv & !sex){
                out<- mFit_var.distri(fecundity=b, j=j, a=a, varJ=varJ, N0=N0)
              }
              if (!broodMortality & varEnv  & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_var.distri")
                # out<- mFitSex_var.distri(fecundity=b, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              if (!broodMortality & varEnv  & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mFit_varseason.distri")
              }
              if (!broodMortality & varEnv  & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_varseason.distri")
              }
              
              if (broodMortality  & !varEnv & !seasonalEnv & !sex){
                out<- mSurvBV.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0)
              }
              if (broodMortality  & !varEnv & !seasonalEnv &  sex){
                out<- mSurvBVSex.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              if (broodMortality  & !varEnv & seasonalEnv  & !sex){
                # TODO: fix bug when breedFail == 1. result is filled with NA when called from run(model)
                # cat("mSurvBV.distri(broods=", broods, ", b=", b, ", j=", j, ", a=", a, ", breedFail=", breedFail, ", N0=", N0, ")\n")
                
                # tmp<- mSurvBV*?_season?*.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0)
                # print(tmp)
                # logP("a") # Error to debug 
                # return(tmp)
                
                ## TODO 
                warning("## TODO: out<- mSurvBV_season.distri")
                # mSurvBV_season.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0)
              }
              if (broodMortality  & !varEnv & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_season.distri")
                # out<- mSurvBVSex_season.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              if (broodMortality  & varEnv  & !seasonalEnv & !sex){
                out<- mSurvBV_var.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0)
              }
              if (broodMortality  & varEnv  & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_var.distri")
                # out<- mSurvBVSex_var.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              if (broodMortality  & varEnv  & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mSurvBV_varseason.distri")
                # out<- mSurvBV_varseason.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0)
              }
              if (broodMortality  & varEnv  & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_varseason.distri")
                # out<- mSurvBVSex_varseason.distri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0)
              }
              
            }else{ # tf > 1
            
              if (!broodMortality & !varEnv & !seasonalEnv & !sex){
                out<- mFit.tdistri(fecundity=b, j=j, a=a, N0=N0, tf=tf)
              }
              if (!broodMortality & !varEnv & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex.tdistri")
                # out<- mFitSex.tdistri(fecundity=b, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
              if (!broodMortality & !varEnv & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mFit_season.tdistri")
              }
              if (!broodMortality & !varEnv & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_season.tdistri")
              }
              if (!broodMortality & varEnv  & !seasonalEnv & !sex){
                ## TODO 
                warning("## TODO: out<- mFit_var.tdistri")
                # out<- mFit_var.tdistri(fecundity=b, j=j, a=a, varJ=varJ, N0=N0, tf=tf)
              }
              if (!broodMortality & varEnv  & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_var.tdistri")
                # out<- mFitSex_var.tdistri(fecundity=b, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
              if (!broodMortality & varEnv  & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mFit_varseason.tdistri")
              }
              if (!broodMortality & varEnv  & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mFitSex_varseason.tdistri")
              }
              
              if (broodMortality  & !varEnv & !seasonalEnv & !sex){
                out<- mSurvBV.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, tf=tf)
              }
              if (broodMortality  & !varEnv & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex.tdistri")
                # out<- mSurvBVSex.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
              if (broodMortality  & !varEnv & seasonalEnv  & !sex){
                # TODO: fix bug when breedFail == 1. result is filled with NA when called from run(model)
                # cat("mSurvBV.tdistri(broods=", broods, ", b=", b, ", j=", j, ", a=", a, ", breedFail=", breedFail, ", N0=", N0, ")\n")
                
                # tmp<- mSurvBV*?_season?*.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, tf=tf)
                # print(tmp)
                # logP("a") # Error to debug 
                # return(tmp)
                
                ## TODO 
                warning("## TODO: out<- mSurvBV_season.tdistri")
                # mSurvBV_season.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, tf=tf)
              }
              if (broodMortality  & !varEnv & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_season.tdistri")
                # out<- mSurvBVSex_season.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
              if (broodMortality  & varEnv  & !seasonalEnv & !sex){
                ## TODO 
                warning("## TODO: out<- mSurvBV_var.tdistri")
                # out<- mSurvBV_var.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, tf=tf)
              }
              if (broodMortality  & varEnv  & !seasonalEnv &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_var.tdistri")
                # out<- mSurvBVSex_var.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
              if (broodMortality  & varEnv  & seasonalEnv  & !sex){
                ## TODO 
                warning("## TODO: out<- mSurvBV_varseason.tdistri")
                # out<- mSurvBV_varseason.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, tf=tf)
              }
              if (broodMortality  & varEnv  & seasonalEnv  &  sex){
                ## TODO 
                warning("## TODO: out<- mSurvBVSex_varseason.tdistri")
                # out<- mSurvBVSex_varseason.tdistri(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, tf=tf)
              }
            }            
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
  N_t1<- distriSum(recruits, survivors)

  return (N_t1)
}

# N_t+1 = B(n=N_t * fecundity, p=Beta(j, varJ)) + B(n=N_t, p=a)
mFit_var.distri<- function(fecundity, j, a, varJ, N0, logP=FALSE){
  betaPars<- fbeta(mean=j, var=varJ)
  # if (all(is.na(betaPars)))
  ## Recruitment probability conditioned to the number of successful broods
  recruits<- distriBetaBinom(N0 * fecundity, shape1=betaPars$shape1, shape2=betaPars$shape2, logP=logP)
  
  ## Adult survival
  survivors<- distriBinom(size=N0, prob=a, logP=logP)
  N_t1<- distriSum(recruits, survivors)
  
  return (N_t1)
}

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
  N_t1<- distriSum(recruits, survivors)
  
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
  N_t1<- distriSum(recruits, survivors)
  
  return (N_t1)
}

##TODO: what if there is differential mortalities between males and females?
## SEX RATIO: reproductive population size depends on the behavior. Monogamous = min(males, females), polygynious = females, polyandrious = males

# Adult mortality + offspring mortality + sex ratio
# N_t+1 = B(n=B(n=N_t * fecundity, p=j), p=sexRatio) + B(n=N_t, p=a)
mFitSex.distri<- function(fecundity, j, a, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], Nf, Nm=Nf, logP=FALSE){
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
mSurvBVSex.distri<- function(broods, b, j, a, breedFail, sexRatio=0.5, matingSystem=c("monogamy", "polygyny", "polyandry")[1], Nf, Nm, logP=FALSE){

}

## Multiple time steps models ----
## Decrease accuracy after 50 timesteps aprox. Use logP(p) reduce omitted probability
mFit.tdistri<- function(fecundity, j, a, N0, tf, logP=TRUE){
  N<- mFit.distri(fecundity, j, a, N0, logP=logP)
  at<- attributes(N)
  for (i in 2:tf){
    cat(i, "/", tf, nrow(N), "p=", sum(logP(N, logP=FALSE)$p), "max x=", N$x[nrow(N)], "\n")
    N<- mFit.distri(fecundity, j, a, N, logP=logP)
    print(utils::head(N$x[N$p %in% c(-Inf, 0)]))
    N<- N[!N$p %in% c(-Inf, 0),]
  }
  attributes(N)$p.omitted<- 1 - sum(logP(N, logP=FALSE)$p)
  attributes(N)$parameters<- at$parameters
  attributes(N)$tf<- tf
  
  return (N)
}

mSurvBV.tdistri<- function(broods, b, j, a, breedFail, N0, tf, logP=TRUE){
  N<- mSurvBV.distri(broods, b, j, a, breedFail, N0, logP=logP)
  at<- attributes(N)
  for (i in 2:tf){
    cat(i, "/", tf, nrow(N), "p=", sum(logP(N, logP=FALSE)$p), "max x=", N$x[nrow(N)], "\n")
    N<- mSurvBV.distri(broods, b, j, a, breedFail, N, logP=logP)
    N<- N[!N$p %in% c(-Inf, 0),]
  }
  attributes(N)$p.omitted<- 1 - sum(N$p)
  attributes(N)$parameters<- at$parameters
  attributes(N)$tf<- tf
  
  return (N)
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
lambda.numericDistri<- function(distri, N0, tf=1){
  distri$x<- (distri$x / N0)^(1/tf)
  return(distri)
}

r.numericDistri<- function(distri, N0, tf=1){
  distri$x<- (distri$x - N0) / N0 / tf
  return(distri)
}
