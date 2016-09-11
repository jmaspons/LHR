#' @include aaa-classes.R
NULL

## Model Constructor ----

#' Model constructor
#'
#' @rdname Model
#' @param lh 
#' @param env 
#' @param sim 
#' @return a \code{Model} object.
#' @examples model<- Model()
#' @export
setGeneric("Model", function(lh=LH(method="LH axes"), env=Env(), sim=Sim(type=match.arg(type)), pars, type=c("discretePopSim", "numericDistri", "ABM", "ssa"), ...) standardGeneric("Model"))

# c("discretePopSim", "numericDistri", "ABM", "ssa")
setMethod("Model",
          signature(lh="ANY", env="ANY", sim="ANY", pars="missing", type="ANY"),
          function(lh=LH(method="LH axes"), env=Env(), sim=Sim(type=match.arg(type)), type=c("discretePopSim", "numericDistri", "ABM", "ssa"), ...){
            
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
              
              if (inherits(sim, "Sim.discretePopSim")){
                model<- new("Model.discretePopSim", scenario, sim=sim, params=parameters)
              }else{
                model<- new("Model.numericDistri", scenario, sim=sim, params=parameters)
              }
            }
            
            return (model)
          }
)

setMethod("Model",
          signature(lh="missing", env="missing", sim="ANY", pars="data.frame", type="ANY"),
          function(sim=Sim(type=match.arg(type)), pars, type=c("discretePopSim", "numericDistri", "ABM", "ssa")){
            modelClass<- gsub("Sim.", "Model.", class(sim))
            new(modelClass, pars, sim=sim)
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
            tmpLH<- data.frame(idLH=rownames(lh), S3Part(lh), interval, stringsAsFactors=FALSE)
            tmpEnv<- data.frame(idEnv=rownames(env), S3Part(env), stringsAsFactors=FALSE)

            scenario<- merge(tmpLH, tmpEnv)
            scenario<- data.frame(idScenario=paste0("lh", scenario$idLH, "_env", scenario$idEnv), scenario, stringsAsFactors=FALSE)
            scenario<- scenario[naturalsort::naturalorder(scenario$idScenario),]
            rownames(scenario)<- scenario$idScenario
            
            ## Brood failure (e.g. nest predation)
            # breedFail is a proportion of juvenile mortality correlated at the brood level
            # survival: P(j) = P(jbr AND jind) = P(jbr) * P(jind)
            # death:    P(-j) = P(-jbr OR P(-jind | jbr))
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
#' @examples res<- run(model, cl=2)
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
  scenario<- split(scenario, scenario$idScenario)
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
  
  
  stats<- lapply(sim, function(x) x$stats)
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)

  simRes<- model@sim
  S3Part(simRes)<- stats

  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@raw<- rawSim
  }
  if (pars$Ntf){
    Ntf<- lapply(sim, function(x) x$Ntf)
    Ntf<- do.call("rbind", Ntf)
    Ntf<- as.data.frame(Ntf)
    rownames(Ntf)<- paste0(Ntf$idScenario, "_N", Ntf$N0)
    # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  return(simRes)
}


runScenario.discretePopSim<- function (scenario, pars, verbose=FALSE){
  if (verbose){
    message(rownames(scenario), "/", nrow(scenario), "\n")
    print(scenario, row.names=FALSE)
  }
  
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
               dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0),
                             stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  stats<- as.data.frame(stats)
  stats$idScenario<- scenario$idScenario
  stats$N0<- pars$N0
  
  if (pars$Ntf){
    Ntf<- matrix(nrow=length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0), 
                               Ntf=c("idScenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf)
    Ntf$idScenario<- scenario$idScenario
    Ntf$N0<- pars$N0
  }
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
  seasonVar<- seasonOptimCal(env=Env(scenario), nSteps=scenario$broods, interval=scenario$interval)[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)

  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    
    pop<- with(scenario, discretePopSim(broods=broods, b=b, j=jind, a=a, breedFail=1 - jbr,
               varJ=ifelse(pars$envVar$j, scenario$var, 0), varBreedFail=ifelse(pars$envVar$breedFail, scenario$var, 0),
               seasonVar=seasonVar,
               sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf, maxN=pars$maxN))
    
    if (is.null(pop) | all(is.na(pop)) | is.list(pop)){ ## TODO: pop<- list(popF, popM) when mating system != NA -> 2 sexes
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- N0
      }
      
      next
    }
    
    stats[n, -c(1:2)]<- as.numeric(summary(pop))
    
    if (pars$raw){
      rawSim[[n]]<- pop
      names(rawSim)[n]<- N0
    }
    if (pars$Ntf){
      pop<- pop[,ncol(pop)]
      pop[is.na(pop)]<- 0
      Ntf[n, -(1:2)]<- sort(pop)
    }
  }
  
  res<- list(stats=stats)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}



run.numericDistri<- function(model, cl=parallel::detectCores()){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)
  pars<- model@sim@params
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }

  parallel::clusterExport(cl=cl, "pars", envir=environment())
  parallel::clusterEvalQ(cl, library(LHR))

  sim<- parallel::parLapply(cl=cl, scenario, runScenario.numericDistri, pars=pars)
  
  # sim<- lapply(scenario, runScenario.numericDistri, pars=pars)
  #   
  # sim<- list()
  # for (i in seq_along(scenario)){
  #   sim[[i]]<- runScenario.numericDistri(scenario=scenario[[i]], pars=pars)
  # }
  
  
  stats<- lapply(sim, function(x) x$stats)
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats)
  stats<- stats[naturalsort::naturalorder(paste0(stats$idScenario, "|", stats$N0)),]
  
  simRes<- model@sim
  S3Part(simRes)<- stats
  
  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@raw<- rawSim
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  return(simRes)
}


runScenario.numericDistri<- function(scenario, pars, verbose=FALSE){
  if (verbose){
    cat(rownames(scenario), "/", nrow(scenario), "\n")
    print(scenario, row.names=FALSE)
  }
  
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
               dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0),
                             stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  stats<- as.data.frame(stats)
  stats$idScenario<- scenario$idScenario
  stats$N0<- pars$N0
  
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
                                                 sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, tf=pars$tf)) #TODO: add , maxN=pars$maxN
    if (is.null(distri) | all(is.na(distri))){
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- N0
      }
      next
    }
    
    ## TODO 2 sexes models: pop<- list(popF, popM)
    distri<- logP(distri, logP=FALSE)
    distri<- cumP(distri)
    selN0<- which(distri$x == N0)
    statsTmp<- c(increase= 1 - distri$cump[selN0], decrease=distri$cump[selN0-1], stable=distri$p[selN0], extinct=distri$p[1])
    
    distriLambda<- lambda(distri, N0=N0, tf=pars$tf)
    distriR<- r(distri, N0=N0, tf=pars$tf)
    distriLambda<- sdistri(distriLambda)
    distriR<- sdistri(distriR)
    
    resTmp<- c(statsTmp, as.numeric(distriR), as.numeric(distriLambda))
    stats[n, -c(1:2)]<- resTmp
    
    if (pars$raw){
      rawSim[[n]]<- distri
      names(rawSim)[n]<- N0
    }
  }
  
  res<- list(stats=stats)
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}


run.ABM<- function(model, cl=parallel::detectCores(), raw, ...){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)

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
  
  sim<- parallel::parLapply(cl=cl, scenario, runScenario.ABM, pars=pars)

  # sim<- lapply(scenario, LHR:::runScenario.ABM, pars=pars)
  # 
  # sim<- list()
  # for (i in seq_along(scenario)){
  #   sim[[i]]<- runScenario.discretePopSim(scenario=scenario[[i]], pars=pars)
  # }
  
  
  stats<- lapply(sim, function(x) x$stats)
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats)

  simRes<- model@sim
  S3Part(simRes)<- stats
  

  if (pars$raw){
    simRes@raw<- lapply(sim, function(x) x$raw)
  }
  if (pars$discretePopSim){
    simRes@discretePopSim<- lapply(sim, function(x) x$discretePopSim)
  }
  if (pars$Ntf){
    Ntf<- lapply(sim, function(x) x$Ntf)
    Ntf<- do.call("rbind", Ntf)
    Ntf<- as.data.frame(Ntf)
    rownames(Ntf)<- paste0(Ntf$idScenario, "_N", Ntf$N0)
    # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }
  
  if (numericCL) parallel::stopCluster(cl)
  
  return (simRes)
}


runScenario.ABM<- function (scenario, pars, verbose=FALSE){
  if (verbose){
    message(rownames(scenario), "/", nrow(scenario), "\n")
    print(scenario, row.names=FALSE)
  }
  
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
               dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", sapply(pars$N0, sum)),
                             stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  stats<- as.data.frame(stats)
  stats$idScenario<- scenario$idScenario
  stats$N0<- sapply(pars$N0, sum)
  
  if (pars$Ntf){
    Ntf<- matrix(NA_real_, nrow=length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", sapply(pars$N0, sum)), 
                               Ntf=c("idScenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf)
    Ntf$idScenario<- scenario$idScenario
    Ntf$N0<- sapply(pars$N0, sum)
  }
  
  if (pars$discretePopSim) pop<- list()
  
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
#   seasonVar<- seasonOptimCal(env=Env(scenario))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
#   if (any(seasonVar !=1)){
#     jindSeason<- scenario$jind * seasonVar
#     jbrSeason<- scenario$jbr * seasonVar
#   }else{
#     jindSeason<- scenario$jind
#     jbrSeason<- scenario$jbr
#   }
  
  for (n in 1:length(pars$N0)){
    N0<- pars$N0[[n]]
    popABM<- discreteABMSim(N0=N0, params=scenario, transitionsFunc=pars$transitionsFunc, replicates=pars$replicates, tf=pars$tf, maxN=pars$maxN, Ntf=!pars$raw & !pars$discretePopSim)
    
    N0<- sum(N0)
    
    if (is.null(popABM) | all(is.na(popABM))){
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- paste0("N", N0)
      }
      
      next
    }
    
    
    stats[n, -c(1:2)]<- summary(popABM) # First columns: idScenario and N0
    
    if (pars$raw){
      rawSim[[n]]<- popABM
      names(rawSim)[n]<- N0
    }
    if (pars$discretePopSim){
      pop[[n]]<- discreteABMSim2discretePopSim(popABM)
      names(pop)[n]<- N0
    }
    if (pars$Ntf){
      popABM<- popABM[,,dim(popABM)[3]]
      popABM<- rowSums(popABM)
      popABM[is.na(popABM)]<- 0
      Ntf[n, -c(1:2)]<- sort(popABM) # First columns: idScenario and N0
    }
  }
  
  res<- list(stats=stats)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
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

  res<- new("Sim.ssa", res$stats, Ntf=res$Ntf, params=model@sim@params, raw=res, N0_Pest=model@sim@N0_Pest)
  # simRes<- model@sim
  # S3Part(simRes)<- res$stats
  # 
  # # if (pars$raw){
  # #   simRes@raw<- lapply(res, function(x) x$raw)
  # # }
  # if (discretePop){
  #   simRes@discretePopSim<- lapply(res$pop, function(x) x$pop)
  # }
  # if (finalPop){
  #   Ntf<- lapply(res, function(x) x$Ntf)
  #   Ntf<- do.call("rbind", Ntf)
  #   Ntf<- as.data.frame(Ntf)
  #   rownames(Ntf)<- paste0(Ntf[,"idScenario"], "_N", Ntf[,"N0"])
  #   # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
  #   simRes@Ntf<- Ntf
  # }
  # 
  
  if (numericCL) parallel::stopCluster(cl)
  
  return (res)
}




## result(): Post process result ----
#' @rdname Model
#'
#' @param model 
#' @param type 
#' @details TODO: type="Ntf" doesn't work for Model.ssa. Check it and standardize model@sim@Ntf
#'
#' @return a data frame with the aggregated results and parameters of a simulation.
#' @examples result(res) 
#' @export
setGeneric("result", function(model, type=c("stats", "N0_Pest", "Ntf")) standardGeneric("result"))

setMethod("result", 
          signature(model="Model", type="ANY"),
          function(model, type=c("stats", "N0_Pest", "Ntf")){
            type<- match.arg(type)
            
            if (nrow(model@sim) == 0 & nrow(model@sim@N0_Pest) == 0){
              stop("There are no results yet. Use run(model) to start the simulations. The function return a model with the results on the sim slot.\n")
            }

            res<- switch(type,
              stats={
                if (nrow(model@sim) == 0){
                  stop("There are no results yet. Use run(model) to start the simulations. The function return a model with the results on the sim slot.\n")
                }
                res<- merge(S3Part(model), S3Part(model@sim), by="idScenario")

                # sort columns idScenario, N0, ...
                res<- cbind(res[,c("idScenario","N0")], res[,-c(1, grep("N0", names(res)))]) 
                res<- res[naturalsort::naturalorder(paste0(res$idScenario, "|", res$N0)),]
                res
              },
              N0_Pest={
                if (nrow(model@sim@N0_Pest) == 0){
                  stop("There are no results yet. Use findN0_Pest(model) to start the simulations. The function return a model with the results on the model@sim@N0_Pest slot.\n")
                }
                res<- cbind(S3Part(model), S3Part(model@sim@N0_Pest))
              },
              Ntf={
                if (nrow(model@sim@Ntf) == 0){
                  stop("There are no results yet. Use run(model) to start the simulations. Check that model@sim@params@Ntf is TRUE before running. The function return a model with the results on the @sim@Ntf slot.\n")
                }
                res<- model@sim@Ntf
                res<- reshape2::melt(res, id.vars=1:2, value.name="Ntf") # id vars: idScenario and N0
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
setMethod("show", signature(object="Model"),
          function(object){
            cat("Object of class ", class(object), " with", nrow(object), "scenarios\n\n")
            
            if (nrow(object@sim) == 0 & nrow(object@sim@N0_Pest) == 0){
              print(S3Part(object)) # S3Part(x)[,1:12]
              cat("There are no results yet. Use run(model) to start the simulations for all scenarios and N0 or findN0_Pest(model) to find the N0 giving a fixed extinction probability.\n")
            } else {
            
              if(nrow(object@sim) > 0) {
                cat("Stats:\n")
                res<- result(object)
                print(res)
                cat("Use result(model) to get a data.frame with the parameters and the results.\n\n")
              }
              
              if (nrow(object@sim@N0_Pest) > 0){
                cat("N0_Pest:\n")
                res<- result(object, type="N0_Pest")
                print(res)
                cat("Use result(model, type=\"N0_Pest\") to get a data.frame with the parameters and the results.\n\n")
              }
              
              if (nrow(object@sim) == 0){
                cat("Use run(model) to start the simulations for all scenarios and N0.\n")
              }
              
              if (nrow(object@sim@N0_Pest) == 0){
                cat("Use findN0_Pest(model) to find the N0 giving a fixed extinction probability (0.5 by default).\n")
              }
            }
            
            invisible(object)
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
# TODO: keep the results from run(model)
#' @rdname Model
#' @export
`[.Model`<- function(x, ...){
  xSel<- data.frame(x)[...]
  Model(lh=LH(xSel), env=Env(xSel), sim=Sim(x))
}

#' Plot
#'
#' @rdname Model
#' @param x 
#' @param resultType 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.Model<- function(x, resultType=c("stats", "N0_Pest"), ...){
  resultType<- match.arg(resultType)
  
  if (nrow(x@sim) == 0 & nrow(x@sim@N0_Pest) == 0){
    message("No results found. Plotting the parameter space.\n\trun(model) to simulate.")
    selNum<- sapply(S3Part(x), is.numeric)
    return(graphics::plot(S3Part(x)[, selNum]))
  }
  
  if (nrow(x@sim) == 0 | resultType == "N0_Pest"){  
    res<- result(x, type="N0_Pest")
    res$idScenario<- factor(res$idScenario)
    
    out<- ggplot2::ggplot(res, ggplot2::aes(x=lambda, y=(N0interpoled), group=idScenario, color=idScenario)) + 
      ggplot2::geom_point() + ggplot2::facet_grid(breedFail~seasonAmplitude + var, labeller=ggplot2::label_both)  
    return(out)
  }
  
  if (nrow(x@sim@N0_Pest) == 0 | resultType == "stats")
  res<- result(x, type="stats")
  res$Pest<- 1 - res$extinct
  res$idScenario<- factor(res$idScenario)
  
  out<- ggplot2::ggplot(res, ggplot2::aes(x=N0, y=Pest, group=idScenario, color=idScenario)) + 
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(breedFail~seasonAmplitude + var, labeller=ggplot2::label_both)
  return(out)
}


#' Histogram
#'
#' @rdname Model
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
hist.Model<- function(x, ...){
  if (nrow(x@sim) == 0){
    stop("No results found.\n\trun(model) to simulate.")
    selNum<- sapply(S3Part(x), is.numeric)
    return(graphics::hist(S3Part(x)[, selNum]))
  }
  
  N<- x@sim@Ntf
  xd<- data.frame(x)
  
  N<- merge(N, xd, by="idScenario")
  N$idScenario<- factor(res$idScenario)
  N<- reshape2::melt(N, id.vars=c("idScenario", "N0", "seasonAmplitude", "var", "breedFail"), value.name="Ntf")
  
  ggplot2::ggplot(N, ggplot2::aes(x=Ntf, group=idScenario, color=idScenario)) + 
    ggplot2::geom_histogram() + ggplot2::facet_grid(N0 + breedFail~seasonAmplitude + var, labeller=ggplot2::label_both)
}
