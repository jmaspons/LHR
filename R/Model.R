#' @include aaa-classes.R
NULL

## Model Constructor ----

#' Model constructor
#'
#' @rdname Model
#' @param lh 
#' @param env 
#' @param sim 
#'
#' @return a \code{Model} object.
#' @examples model<- Model()
#' @export
setGeneric("Model", function(lh=LH(), env=Env(), sim=Sim(), pars, ...) standardGeneric("Model"))

setMethod("Model",
          signature(lh="ANY", env="ANY", sim="ANY", pars="missing"),
          function(lh=LH(), env=Env(), sim=Sim(), ...){
            if (inherits(sim, "Sim.ABM")){
              pars<- getParamsCombination.LH_Beh(lh=lh, env=env, ...)
              model<- new("Model.ABM", pars, sim=sim)
            }else if (inherits(sim, "Sim.ssa")){
              pars<- getParamsCombination.LH_Beh.ssa(...)
              model<- new("Model.ssa", pars, sim=sim)
            } else if (inherits(sim, c("Sim.discretePopSim", "Sim.numericDistri"))){
              ## WARNING: Sim.ABM and Sim.ssa inherits from Sim.discretePopSim. Exclude from here
              lhEnv<- combineLH_Env(lh=lh, env=env)
              scenario<- lhEnv$scenario
              parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv) #, breedFail=lhEnv$breedFail)
              
              model<- new("Model", scenario, sim=sim, params=parameters)
              
            }
            
            return (model)
          }
)


setMethod("Model",
          signature(lh="missing", env="missing", sim="Sim.ABM", pars="data.frame"),
          function(sim, pars){
            new("Model.ABM", pars, sim=sim)
          }
)

setMethod("Model",
          signature(lh="missing", env="missing", sim="Sim.ssa", pars="data.frame"),
          function(sim, pars){
            new("Model.ssa", pars, sim=sim)
          }
)


## Combine LH and Environment ----
#' Combine a LH and a Env object in a set of scenarios
#'
#' @param lh 
#' @param env 
#' @param resolution the number of divisions of the sinusoidal pattern representing one year. Used in \code{\link{seasonOptimCal}}.
#' @param interval a vector with the number of events for each LH strategy.
#' @param criterion 
#'
#' @return the nSteps values separed by interval units of a sinusoidal pattern optimizing 
#' according to different criterion (maxMean, maxFirst)
#' seasonEvents<- function(env)
#' @examples combineLH_Env()
#' @export
setGeneric("combineLH_Env", function(lh=LH(), env=Env(), resolution=12, interval=2, criterion=c("maxFirst", "maxMean")[1]) standardGeneric("combineLH_Env"))

setMethod("combineLH_Env", 
          signature(lh="ANY", env="ANY", resolution="ANY", interval="ANY", criterion="ANY"),
          function(lh=LH(), env=Env(), resolution=12, interval=2, criterion="maxFirst"){
            tmpLH<- data.frame(S3Part(lh), interval)
            scenario<- merge(tmpLH, S3Part(env))
            
            ## Brood failure (e.g. nest predation)
            # breedFail is a proportion of juvenile mortality correlated at the brood level
            # survival: P(j) = P(jbr ∩ jind) = P(jbr) * P(jind)
            # death:    P(-j) = P(-jbr ∪ P(-jind | jbr))
            scenario$jind<- scenario$j / (scenario$breedFail * (scenario$j-1) + 1)     # juvenile mortality is divided in brood mortality and individual mortality
            scenario$jbr<- scenario$breedFail * (scenario$j-1) + 1

            ## Seasonality
            seasonBroods<- data.frame(scenario[,c("seasonMean", "seasonAmplitude", "broods", "interval")], var=0, breedFail=0) # var necessary for Env(seasonBroods)
            seasonVar<- seasonOptimCal(env=seasonBroods, resolution=resolution, nSteps=seasonBroods$broods, interval=seasonBroods$interval, criterion=criterion)
            
            return (list(scenario=scenario, seasonBroodEnv=list(parsBroodEnv=seasonBroods, broodEnv=seasonVar))) #, breedFail=breedFail))
          }
)


## run(): Simulate models ----

#' @rdname Model
#'
#' @param model 
#' @param cl 
#' @param ... 
#'
#' @return a \code{Model} object with the result of the simulation.
#' @examples res<- run(model)
#' @export
setGeneric("run", function(model, cl=parallel::detectCores(), ...) standardGeneric("run"))

setMethod("run", 
          signature(model="Model", cl="ANY"),
          function(model, cl=parallel::detectCores(), ...){
            if (is.numeric(cl)){
              numericCL<- TRUE
              cl<- parallel::makeCluster(cl)
            } else {
              numericCL<- FALSE
            }
            
            simRes<- switch(class(model@sim),
                               Sim.discretePopSim=run.discretePopSim(model, cl=cl),
                               Sim.numericDistri=run.numericDistri(model, cl=cl),
                               Sim.ABM=run.ABM(model, cl=cl, ...),
                               Sim.ssa=run.ssa(model, cl=cl, ...))

            modelRes<- new("Model", 
                        S3Part(model),
                        sim=simRes,
                        params=model@params)
            
            if (numericCL) parallel::stopCluster(cl)
            
            return(modelRes)
          }
          
)


run.discretePopSim<- function(model, cl=parallel::detectCores()){
  scenario<- S3Part(model)
  scenario<- split(scenario, rownames(scenario))
  pars<- model@sim@params

  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }
  
  parallel::clusterExport(cl=cl, "pars", envir=environment())
  parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
  parallel::clusterEvalQ(cl, library(LHR))
  sim<- parallel::parLapply(cl=cl, scenario, runScenario.discretePopSim, pars=pars)

  # sim<- lapply(scenario, LHR:::runScenario.discretePopSim, pars=pars)
  # 
  # sim<- list()
  # for (i in seq_along(scenario)){
  #   sim[[i]]<- runScenario.discretePopSim(scenario=scenario[[i]], pars=pars)
  # }
  
  
  res<- lapply(sim, function(x) x$res)
  res<- do.call("rbind", res)
  res<- as.data.frame(res)

  simRes<- new("Sim.discretePopSim", res, params=pars)

  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@raw<- rawSim
  }
  if (pars$Ntf){
    Ntf<- lapply(sim, function(x) x$Ntf)
    Ntf<- do.call("rbind", Ntf)
    Ntf<- as.data.frame(Ntf)
    rownames(Ntf)<- paste0(Ntf$scenario, "_N", Ntf$N0)
    # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  return(simRes)
}


runScenario.discretePopSim<- function (scenario, pars){
  message(rownames(scenario), "/", nrow(scenario), "\n")
  print(scenario, row.names=FALSE)
  
  res<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
               dimnames=list(scenario=paste0(rep(rownames(scenario), length=length(pars$N0)), "_N", pars$N0),
                             stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  if (pars$Ntf){
    Ntf<- matrix(nrow=length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario=paste0(rep(rownames(scenario), length=length(pars$N0)), "_N", pars$N0), 
                               Ntf=c("scenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf)
  }
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
  seasonVar<- seasonOptimCal(env=Env(scenario))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
  if (any(seasonVar !=1)){
    jindSeason<- scenario$jind * seasonVar
    jbrSeason<- scenario$jbr * seasonVar
  }else{
    jindSeason<- scenario$jind
    jbrSeason<- scenario$jbr
  }
  
  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    
    pop<- with(scenario, discretePopSim(broods=broods, b=b, j=jindSeason, a=a, breedFail=1 - jbrSeason,
               varJ=ifelse(pars$envVar$j, var, 0), varBreedFail=ifelse(pars$envVar$breedFail, var, 0),
               sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf))
    
    if (is.null(pop) | is.na(pop) | is.list(pop)){ ## TODO: pop<- list(popF, popM) when mating system != NA -> 2 sexes
      res[n, c("scenario", "N0")]<- c(as.numeric(rownames(scenario)), N0)
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- N0
      }
      if (pars$Ntf){
        Ntf[n, c("scenario", "N0")]<- c(as.numeric(rownames(scenario)), N0)
      }
      next
    }
    
    res[n,]<- c(as.numeric(rownames(scenario)), N0, as.numeric(summary(pop)))
    
    if (pars$raw){
      rawSim[[n]]<- pop
      names(rawSim)[n]<- N0
    }
    if (pars$Ntf){
      pop<- pop[,ncol(pop)]
      pop[is.na(pop)]<- 0
      Ntf[n,]<- c(as.numeric(rownames(scenario)), N0, sort(pop))
    }
  }
  
  res<- list(res=res)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}



run.numericDistri<- function(model, cl=parallel::detectCores()){
  scenario<- S3Part(model)
  scenario<- split(scenario, rownames(scenario))
  pars<- model@sim@params
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }

  parallel::clusterExport(cl=cl, "pars", envir=environment())
  parallel::clusterEvalQ(cl, library(LHR))
  sim<- parallel::parLapply(cl=cl, X=scenario, fun=runScenario.numericDistri, pars=pars)
  
#   sim<- lapply(scenario, runScenario.numericDistri, pars=pars)
#   
#   sim<- list()
#   for (i in seq_along(scenario)){
#     sim[[i]]<- runScenario.numericDistri(scenario=scenario[[i]], pars=pars)
#   }
  
  
  res<- lapply(sim, function(x) x$res)
  res<- do.call("rbind", res)
  res<- as.data.frame(res)
  res<- res[order(res$scenario, res$N0),]
  
  simRes<- new("Sim.discretePopSim", res, params=pars)
  
  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@raw<- rawSim
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  return(simRes)
}


runScenario.numericDistri<- function(scenario, pars){
  cat(rownames(scenario), "/", nrow(scenario), "\n")
  print(scenario, row.names=FALSE)
  
  res<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
               dimnames=list(scenario=paste0(rep(rownames(scenario), length=length(pars$N0)), "_N", pars$N0),
                             stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
  seasonVar<- seasonOptimCal(env=Env(scenario))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
  if (any(seasonVar !=1)){
    jindSeason<- scenario$jind * seasonVar
    jbrSeason<- scenario$jbr * seasonVar
  }else{
    jindSeason<- scenario$jind
    jbrSeason<- scenario$jbr
  }
  
  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    ## TODO: check error when breedFail == 1:
    # only happens on run(model), not when calling mSurvBV.distri with the same parameters!!
    # maybe var collide with var()?? It seems is not the case...
    distri<- with(scenario, tDistri(broods=broods, b=b, j=jindSeason, a=a, breedFail=1 - jbrSeason,
                                                 varJ=ifelse(pars$envVar$j, var, 0), varBreedFail=ifelse(pars$envVar$breedFail, var, 0),
                                                 sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, tf=pars$tf))
    if (is.null(distri) | is.na(distri)){
      res[n,c("scenario", "N0")]<- c(as.numeric(rownames(scenario)), N0)
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- N0
      }
      next
    }
    
    ## TODO: pop<- list(popF, popM)
    distri<- logP(distri, logP=FALSE)
    distri<- cumsum(distri)
    selN0<- which(distri$x == N0)
    tmp<- c(increase= 1 - distri$cump[selN0], decrease=distri$cump[selN0-1], stable=distri$p[selN0], extinct=distri$p[1])
    distriLambda<- lambda(distri, N0=N0, tf=pars$tf)
    distriR<- r(distri, N0=N0, tf=pars$tf)
    
    res[n,]<- c(as.numeric(rownames(scenario)), N0, tmp, as.numeric(sdistri(distriR)), as.numeric(sdistri(distriLambda)))
    
    if (pars$raw){
      rawSim[[n]]<- distri
      names(rawSim)[n]<- N0
    }
  }
  
  res<- list(res=res)
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}


run.ABM<- function(model, cl=parallel::detectCores(), raw, ...){
  x0L<- model@sim@params$N0
  params<- S3Part(model)
  transitionsFunc<- model@sim@params$transitionsFunc
  tf<- model@sim@params$tf
  replicates<- model@sim@params$replicates
  raw<- model@sim@params$raw
  discretePop<- model@sim@params$discretePopSim
  finalPop<- model@sim@params$Ntf
  #   burnin=-1
  #   dtDiscretize=NULL
  #   cl=1
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }
  
  res<- exploreABM(x0L=x0L, params=params, transitionsFunc=transitionsFunc, 
                   tf=tf, replicates=replicates, raw=raw, discretePop=discretePop, finalPop=finalPop, cl=cl, ...)
  
  out<- new("Sim.ABM", res$stats, params=model@sim@params)
  
  if (finalPop)    out@Ntf           <- res$Ntf
  if (discretePop) out@discretePopSim<- res$pop
  if (raw)         out@raw           <- res$discreABMSim
  
  if (numericCL) parallel::stopCluster(cl)
  
  return (out)
}


run.ssa<- function(model, cl=parallel::detectCores(), ...){
  x0L<- model@sim@params$N0
  params<- S3Part(model)
  transitionMat<- model@sim@params$transitionMat
  rateFunc<- model@sim@params$rateFunc
  tf<- model@sim@params$tf
  replicates<- model@sim@params$replicates
  discretePop<- model@sim@params$raw
  finalPop<- model@sim@params$Ntf
  #   burnin=-1
  #   dtDiscretize=NULL
  #   cl=1
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }
  
  res<- exploreSSA(x0L=x0L, params=params, transitionMat=transitionMat, rateFunc=rateFunc, 
                   maxTf=tf, replicates=replicates, discretePop=discretePop, finalPop=finalPop, cl=cl, ...)

  res<- new("Sim.ssa", res$stats, Ntf=res$Ntf, params=model@sim@params, raw=res)
  
  if (numericCL) parallel::stopCluster(cl)
  
  return (res)
}




## Post process result ----
#' @rdname Model
#'
#' @param model 
#' @param type 
#' @details TODO: type="Ntf" doesn't work for Model.ssa. Check it and standardize model@sim@Ntf
#'
#' @return a data frame with the aggregated results and parameters of a simulation.
#' @examples result(res) 
#' @export
setGeneric("result", function(model, type="stats") standardGeneric("result"))

setMethod("result", 
          signature(model="Model", type="ANY"),
          function(model, type="stats"){
            if (nrow(model@sim) == 0){
              stop("There are no results yet. Use run(model) to start the simulations. The function return a model with the results on the sim slot.\n")
            }

            res<- switch(type,
              stats={
                res<- merge(S3Part(model), S3Part(model@sim), by.x="row.names", by.y="scenario")
                names(res)[1]<- "scenario"
                if (is.numeric(model@sim$scenario)) res$scenario<- as.numeric(res$scenario)
                res<- cbind(res[,c("scenario","N0")], res[,-c(1, grep("N0", names(res)))]) # sort columns scenario, N0, ...
                res<- res[order(res$scenario, res$N0),]
                res},
              Ntf={
                res<- model@sim@Ntf
                res<- reshape2::melt(res, id.vars=1:2, value.name="Ntf") # id vars: scenario and N0
                res$quantile<- as.numeric(res$variable)
                res$quantile<- (res$quantile - 1) / (model@sim@params$replicates - 1) # length(unique(res$quantile)) == replicates
                res<- res[,-3]
              }
            )

          return(res)
        }
)


## Generic methods ----
#' @export
setMethod("print", signature(x="Model"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)

#' @export
setMethod("show", signature(object="Model"),
          function(object){
            cat("Object of class ", class(object), " with", nrow(object), "scenarios\n")
            if (nrow(object@sim) == 0){
              print(S3Part(object)) # S3Part(x)[,1:12]
              cat("There are no results yet. Use run(model) to start the simulations.\n")
            }else{
              res<- result(object)
              print(res)
              cat("Use result(model) to get a data.frame with the parameters and the results.\n")
            }
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
#' @rdname Model
#' @export
`[.Model`<- function(x, ...){
  xSel<- data.frame(x)[...]
  Model(lh=LH(xSel), env=Env(xSel), sim=Sim(x))
}


