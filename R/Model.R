#' @include aaa-classes.R
NULL

## Model Constructor ----

#' Model constructor
#'
#' @rdname Model
#' @param lh a \code{\link{LH}} object. 
#' @param env a \code{\link{Env}} object. 
#' @param sim 
#'
#' @return a \code{Model} object.
#' @examples model<- Model()
#' @export
setGeneric("Model", function(lh=LH(method="LH axes"), env=Env(), sim=Sim(type=match.arg(type)), pars, type=c("discretePopSim", "numericDistri", "ABM", "numericDistriABM"), ...) standardGeneric("Model"))

# c("discretePopSim", "numericDistri", "ABM")
setMethod("Model",
          signature(lh="ANY", env="ANY", sim="ANY", pars="missing", type="ANY"),
          function(lh=LH(method="LH axes"), env=Env(), sim=Sim(type=match.arg(type)), type=c("discretePopSim", "numericDistri", "ABM", "numericDistriABM"), ...){
            
            if (inherits(sim, c("Sim.ABM", "Sim.numericDistriABM"))){
              pars<- getParamsCombination.LHEnv_2patchBeh(lh=lh, env=env, ...)
              if (inherits(sim, "Sim.ABM")){
                model<- new("Model.ABM", pars, sim=sim)
              }else{
                model<- new("Model.numericDistriABM", pars, sim=sim)
              }
            } else if (inherits(sim, c("Sim.discretePopSim", "Sim.numericDistri"))){
              ## WARNING: Sim.ABM inherits from Sim.discretePopSim. Exclude from here
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
          function(sim=Sim(type=match.arg(type)), pars, type=c("discretePopSim", "numericDistri", "ABM")){
            modelClass<- gsub("Sim\\.", "Model.", class(sim))
            type<- gsub("Sim\\.", "", class(sim))
            new(modelClass, pars, sim=sim)
          }
)


## Combine LH and Environment ----
#' Combine a LH and a Env object in a set of scenarios
#'
#' @param lh a \code{\link{LH}} object. 
#' @param env a \code{\link{Env}} object. 
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
            tmpLH<- data.frame(S3Part(lh), interval, stringsAsFactors=FALSE)
            tmpEnv<- S3Part(env)

            scenario<- merge(tmpLH, tmpEnv)
            scenario<- data.frame(idScenario=paste0("LH", scenario$idLH, "_Env", scenario$idEnv), scenario, stringsAsFactors=FALSE)
            
            if ("baseLH" %in% names(scenario)){
              scenario<- scenario[order(scenario$baseLH, scenario$lambda, scenario$idEnv),]
            } else {
              scenario<- scenario[order(scenario$idLH, scenario$idEnv),]
            }
            
            rownames(scenario)<- scenario$idScenario
            
            ## Brood failure (e.g. nest predation)
            # breedFail is a proportion of juvenile mortality correlated at the brood level
            # survival: P(j) = P(jbr AND jind) = P(jbr) * P(jind)
            # death:    P(-j) = P(-jbr OR P(-jind | jbr))
            scenario$jind<- scenario$j / (scenario$breedFail * (scenario$j-1) + 1)     # juvenile mortality is divided in brood mortality and individual mortality
            scenario$jbr<- scenario$breedFail * (scenario$j-1) + 1
            
            ## Check combinations of mean and variance falling outside the domain of the Beta distribution
            betaParsJind<- fbeta(mean=scenario$jind, var=scenario$varJ)
            betaParsJbr<- suppressWarnings(fbeta(mean=scenario$jbr, var=scenario$varJ)) # breedFail == 0 -> jbr == 1: error
            betaParsA<- fbeta(mean=scenario$a, var=scenario$varA)
            
            if (any(errJind<- !is.finite(betaParsJind$shape1) & scenario$varJ != 0) |
                any(errJbr<- !is.finite(betaParsJbr$shape1)  & scenario$varJ != 0 & scenario$breedFail > 0) |
                any(errA<- !is.finite(betaParsA$shape1)    & scenario$varA != 0) ){
              err<- errJind | errJbr | errA
              scenario$errorBeta<- err
              
              if (any(errJind)){
                scenario$maxVarJind<- NA
                betaPars<- fbeta(mean=scenario$jind[errJind], var="max")
                scenario$maxVarJind[errJind]<- sbeta(betaPars$shape1, betaPars$shape2)$var
              }
              if (any(errJbr)){
                scenario$maxVarJbr<- NA
                betaPars<- fbeta(mean=scenario$jbr[errJbr], var="max")
                scenario$maxVarJbr[errJbr]<- sbeta(betaPars$shape1, betaPars$shape2)$var
              }
              if (any(errA)){
                scenario$maxVarA<- NA
                betaPars<- fbeta(mean=scenario$a[errA], var="max")
                scenario$maxVarA[errA]<- sbeta(betaPars$shape1, betaPars$shape2)$var
              }
              
              errJ<- any(errJind | errJbr) 
              errA<- any(errA)
              w<- ""
              if (errJ) w<- "varJ"
              if (errA) w<- paste0(w, "varA", sep=" and ")
              warning("Some combinations of ",  w, " with the mean survival fall outside the parameter space of the Beta distribution. Check errorBeta and maxVar column.")
            }

            ## Seasonality
            seasonBroods<- data.frame(scenario[,c("seasonMean", "seasonAmplitude", "broods", "interval")], varJ=0, varA=0, breedFail=0) # varJ & varA necessary for Env(seasonBroods)
            seasonVar<- seasonOptimCal(env=seasonBroods, resolution=resolution, nSteps=seasonBroods$broods, interval=seasonBroods$interval, criterion=criterion)
            
            return (list(scenario=scenario, seasonBroodEnv=list(parsBroodEnv=seasonBroods, broodEnv=seasonVar))) #, breedFail=breedFail))
          }
)


## run(): Simulate models ----

#' @rdname Model
#'
#' @param model 
#' @param cl The number of cores to use or a cluster object (\code{\link[parallel]{makeCluster}} or 
#'   \code{\link[snow]{makeCluster} from \href{https://cran.r-project.org/web/packages/snow/index.html}{snow} package}) 
#' @param pb if \code{TRUE} and \link[pbapply]{pbapply} package is installed, show a progress bar.
#' @param debug if \code{TRUE} run the simulations in a simple loop and print information about the state.
#' @param ... 
#'
#' @return a \code{Model} object with the result of the simulation.
#' @examples res<- run(model, cl=2)
#' @export
setGeneric("run", function(model, cl=parallel::detectCores(), pb=FALSE, debug=FALSE,...) standardGeneric("run"))

setMethod("run", 
          signature(model="Model", cl="ANY", pb="ANY"),
          function(model, cl=parallel::detectCores(), pb=FALSE, ...){
            if (is.numeric(cl)){
              if (.Platform$OS.type == "windows"){
                cl<- parallel::makePSOCKcluster(cl)
              }else{
                cl<- parallel::makeForkCluster(cl)
              }
              on.exit(parallel::stopCluster(cl))
            }
            
            simRes<- switch(class(model@sim),
                               Sim.discretePopSim=run.discretePopSim(model, cl=cl, pb=pb, debug=debug, ...),
                               Sim.numericDistri=run.numericDistri(model, cl=cl, pb=pb, debug=debug, ...),
                               Sim.ABM=run.ABM(model, cl=cl, pb=pb, debug=debug, ...),
                               Sim.numericDistriABM=run.numericDistriABM(model, cl=cl, pb=pb, debug=debug, ...))

            modelRes<- new(class(model), 
                        S3Part(model),
                        sim=simRes,
                        params=model@params)
            
            return(modelRes)
          }
          
)


run.discretePopSim<- function(model, cl=parallel::detectCores(), pb=FALSE, debug=FALSE){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)
  pars<- model@sim@params
  
  if (!debug){
    if (is.numeric(cl)){
      cl<- parallel::makeCluster(cl)
      on.exit(parallel::stopCluster(cl))
    }
    
    parallel::clusterExport(cl=cl, "pars", envir=environment())
    parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
    parallel::clusterEvalQ(cl, library(LHR))
    
    if (pb & requireNamespace("pbapply", quietly=TRUE)){
      sim<- pbapply::pblapply(scenario, runScenario.discretePopSim, pars=pars, cl=cl)
    }else{
      sim<- parallel::parLapply(cl=cl, scenario, runScenario.discretePopSim, pars=pars)
    }
    
  } else {
    sim<- list()
    for (i in seq_along(scenario)){
      message(i, "/", length(scenario))
      print(scenario[[i]], row.names=FALSE)
      eTime<- system.time(sim[[i]]<- runScenario.discretePopSim(scenario=scenario[[i]], pars=pars))
      message("Elapsed time: ", eTime["elapsed"])
    }
  }
  
  stats<- lapply(sim, function(x) x$stats)
  names(stats)<- NULL
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)

  simRes<- model@sim
  S3Part(simRes)<- stats

  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@discretePopSim<- rawSim
  }
  if (pars$Ntf){
    Ntf<- lapply(sim, function(x) x$Ntf)
    Ntf<- do.call("rbind", Ntf)
    Ntf<- as.data.frame(Ntf, stringsAsFactors=FALSE)
    rownames(Ntf)<- paste0(Ntf$idScenario, "_N", Ntf$N0)
    # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }
  
  return(simRes)
}


runScenario.discretePopSim<- function (scenario, pars){
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=15, 
               dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0),
                             stats=c("idScenario", "N0", "increase", "decrease", "stable","extinct", 
                             "increaseTrans", "decreaseTrans", "stableTrans", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 15 = ncol(summary(pop)) + 2 (id,
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)
  stats$idScenario<- scenario$idScenario
  stats$N0<- pars$N0
  
  if (pars$Ntf){
    Ntf<- matrix(NA_real_, nrow=length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0), 
                               Ntf=c("idScenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf, stringsAsFactors=FALSE)
    Ntf$idScenario<- scenario$idScenario
    Ntf$N0<- pars$N0
  }
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
  seasonVar<- seasonOptimCal(env=Env(scenario), nSteps=scenario$broods, interval=scenario$interval)[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)

  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    
    pop<- with(scenario, discretePopSim(broods=broods, b=b, j=jind, s=s, a=a, AFR=AFR, breedFail=1 - jbr,
               varJ=ifelse(pars$envVar$j, scenario$varJ, 0), varBreedFail=ifelse(pars$envVar$breedFail, scenario$varJ, 0),
               varA=varA, seasonVar=seasonVar,
               sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf, maxN=pars$maxN))
    
    if (is.null(pop) | all(is.na(pop)) | is.list(pop)){ ## TODO: pop<- list(popF, popM) when mating system != NA -> 2 sexes
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- N0
      }
      
      next
    }
    
    stats[n, -c(1:2)]<- summary(pop)
    
    if (pars$raw){
      rawSim[[n]]<- pop
      names(rawSim)[n]<- N0
    }
    if (pars$Ntf){
      tmpNtf<- apply(pop[, -1, drop=FALSE], 1, function(x) x[which(match(x, NA) == 1)[1] - 1])
      tmpNtf[is.na(tmpNtf)]<- pop[is.na(tmpNtf), ncol(pop)] # replicates which run until tf (no extinction nor maxN)
      Ntf[n, -(1:2)]<- sort(tmpNtf)
    }
  }
  
  res<- list(stats=stats)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}


run.numericDistri<- function(model, cl=parallel::detectCores(), pb=FALSE, debug=FALSE){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)
  pars<- model@sim@params
  
  if (!debug){
    if (is.numeric(cl)){
      cl<- parallel::makeCluster(cl)
      on.exit(parallel::stopCluster(cl))
    }
  
    parallel::clusterExport(cl=cl, "pars", envir=environment())
    parallel::clusterEvalQ(cl, library(LHR))
  
    if (pb & requireNamespace("pbapply", quietly=TRUE)){
      sim<- pbapply::pblapply(scenario, runScenario.numericDistri, pars=pars, cl=cl)
    }else{
      sim<- parallel::parLapply(cl=cl, scenario, runScenario.numericDistri, pars=pars)
    }
    
  } else {
    sim<- list()
    for (i in seq_along(scenario)){
      message(i, "/", length(scenario))
      print(scenario[[i]], row.names=FALSE)
      eTime<- system.time(sim[[i]]<- runScenario.numericDistri(scenario=scenario[[i]], pars=pars))
      message("Elapsed time: ", eTime["elapsed"])
    }
  }
  
  stats<- lapply(sim, function(x) x$stats)
  names(stats)<- NULL
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)

  simRes<- model@sim
  S3Part(simRes)<- stats
  
  if (pars$raw){
    rawSim<- lapply(sim, function(x) x$raw)
    simRes@raw<- rawSim
  }

  return(simRes)
}


runScenario.numericDistri<- function(scenario, pars){
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=12, 
            dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", pars$N0),
                          stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = lenght(summary(pop)) + 2 (idScenario + N0)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)
  stats$idScenario<- scenario$idScenario
  stats$N0<- pars$N0
  
  if (pars$raw) rawSim<- list()
  
  ## Seasonality
  seasonVar<- seasonOptimCal(env=Env(scenario))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)

  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    ## TODO: check error when breedFail == 1:
    # only happens on run(model), not when calling mSurvBV.distri with the same parameters!!
    # maybe var collide with var()?? It seems is not the case...
    distri<- with(scenario, tDistri(broods=broods, b=b, j=jind, a=a, breedFail=1 - jbr,
                                   varJ=ifelse(pars$envVar$j, scenario$varJ, 0), varBreedFail=ifelse(pars$envVar$breedFail, scenario$varJ, 0),
                                   seasonVar=seasonVar,
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


run.ABM<- function(model, cl=parallel::detectCores(), raw, pb=FALSE, debug=FALSE, ...){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)
  pars<- model@sim@params

  if (!debug){
    if (is.numeric(cl)){
      cl<- parallel::makeCluster(cl)
      on.exit(parallel::stopCluster(cl))
    }
    
    parallel::clusterExport(cl=cl, "pars", envir=environment())
    parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
    parallel::clusterEvalQ(cl, library(LHR))
    
    if (pb & requireNamespace("pbapply", quietly=TRUE)){
      sim<- pbapply::pblapply(scenario, runScenario.ABM, pars=pars, cl=cl)
    }else{
      sim<- parallel::parLapply(cl=cl, scenario, runScenario.ABM, pars=pars)
    }
    
  } else {
    sim<- list()
    for (i in seq_along(scenario)){
      message(i, "/", length(scenario))
      print(scenario[[i]], row.names=FALSE)
      eTime<- system.time(sim[[i]]<- runScenario.ABM(scenario=scenario[[i]], pars=pars))
      message("Elapsed time: ", eTime["elapsed"])
    }
  }
  
  stats<- lapply(sim, function(x) x$stats)
  names(stats)<- NULL
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)

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
    Ntf<- as.data.frame(Ntf, stringsAsFactors=FALSE)
    rownames(Ntf)<- paste0(Ntf$idScenario, "_N", Ntf$N0)
    # Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }

  return (simRes)
}


runScenario.ABM<- function (scenario, pars, randomizeN0=FALSE){
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=15, 
               dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", sapply(pars$N0, sum)),
                             stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct",
                             "increaseTrans", "decreaseTrans", "stableTrans", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 15 = ncol(summary(pop)) + 2 (id, N0)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)
  stats$idScenario<- scenario$idScenario
  stats$N0<- sapply(pars$N0, sum)
  
  if (pars$Ntf){
    Ntf<- matrix(NA_real_, nrow=length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", sapply(pars$N0, sum)), 
                               Ntf=c("idScenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf, stringsAsFactors=FALSE)
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
    popABM<- discreteABMSim(N0=N0, params=scenario, transitionsFunc=pars$transitionsFunc, replicates=pars$replicates,
                            tf=pars$tf, maxN=pars$maxN, Ntf=!pars$raw & !pars$discretePopSim, randomizeN0=randomizeN0)
    
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
      if (pars$discretePopSim) popTmp<- pop[[n]]
      else popTmp<- discreteABMSim2discretePopSim(popABM)
      
      if (ncol(popTmp) == 2){ ## Only one timestep or (pars$Ntf & !pars$raw & !pars$discretePopSim)
        Ntf[n, -(1:2)]<- sort(popTmp[, 2]) # First columns: idScenario and N0
        
      } else {
        tf<- apply(popTmp, 1, function(x) which(match(x, NA) == 1)[1] - 1)
        tf[is.na(tf)]<- ncol(popTmp)
        
        # check if tf is maxN
        if (length(tfMaxN<- unique(tf)) == 1 & all(popTmp[, tfMaxN] > 0)){
          Ntf[n, -(1:2)]<- sort(popTmp[, tfMaxN - 1]) # take last value before saturation (useful for ABM where tfMaxN could be much larger than the average N)
        } else {
          tmpNtf<- mapply(function(x, tf) x[tf], x=split(popTmp, 1:nrow(popTmp)), tf=tf)
          Ntf[n, -(1:2)]<- sort(tmpNtf) # First columns: idScenario and N0
        }
      }
    }
  }
  
  res<- list(stats=stats)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$discretePopSim) res<- c(res, list(discretePopSim=pop))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}


run.numericDistriABM<- function(model, cl=parallel::detectCores(), raw, pb=FALSE, debug=FALSE, ...){
  scenario<- S3Part(model)
  scenario<- split(scenario, scenario$idScenario)
  pars<- model@sim@params
  
  if (!debug){
    if (is.numeric(cl)){
      cl<- parallel::makeCluster(cl)
      on.exit(parallel::stopCluster(cl))
    }
    
    parallel::clusterExport(cl=cl, "pars", envir=environment())
    parallel::clusterEvalQ(cl, library(LHR))
    
    if (pb & requireNamespace("pbapply", quietly=TRUE)){
      sim<- pbapply::pblapply(scenario, runScenario.numericDistriABM, pars=pars, cl=cl)
    }else{
      sim<- parallel::parLapply(cl=cl, scenario, runScenario.numericDistriABM, pars=pars)
    }
    
  } else {
    sim<- list()
    for (i in seq_along(scenario)){
      message(i, "/", length(scenario))
      print(scenario[[i]], row.names=FALSE)
      eTime<- system.time(sim[[i]]<- runScenario.numericDistriABM(scenario=scenario[[i]], pars=pars))
      message("Elapsed time: ", eTime["elapsed"])
    }
  }
  
  stats<- lapply(sim, function(x) x$stats)
  names(stats)<- NULL
  stats<- do.call("rbind", stats)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)
  
  simRes<- model@sim
  S3Part(simRes)<- stats
  
  
  if (pars$raw){
    simRes@raw<- lapply(sim, function(x) x$raw)
    names(simRes@raw)<- sapply(scenario, function(x) x$idScenario)
  }
  if (pars$numericDistriSim){
    simRes@numericDistriSim<- lapply(sim, function(x) x$numericDistriSim)
    names(simRes@numericDistriSim)<- sapply(scenario, function(x) x$idScenario)
  }
  if (pars$Ntf){
    simRes@Ntf<- lapply(sim, function(x) x$Ntf)
    names(simRes@Ntf)<- sapply(scenario, function(x) x$idScenario)
  }
  
  return (simRes)
}


runScenario.numericDistriABM<- function (scenario, pars, randomizeN0=FALSE){
  stats<- matrix(NA_real_, nrow=length(pars$N0), ncol=15, 
                 dimnames=list(scenario_N0=paste0(rep(scenario$idScenario, length=length(pars$N0)), "_N", sapply(pars$N0, sum)),
                               stats=c("idScenario", "N0", "increase", "decrease", "stable", "extinct",
                                       "increaseTrans", "decreaseTrans", "stableTrans", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 15 = ncol(summary(pop)) + 2 (id, N0)
  stats<- as.data.frame(stats, stringsAsFactors=FALSE)
  stats$idScenario<- scenario$idScenario
  stats$N0<- sapply(pars$N0, sum)
  
  if (pars$Ntf){
    Ntf<- list()
  }
  
  if (pars$numericDistriSim) distriSim<- list()
  
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
    distriABM<- numericDistriABMSim(N0=N0, params=scenario, transitionsFunc=pars$transitionsFunc,
                            tf=pars$tf, maxN=pars$maxN, Ntf=!pars$raw & !pars$numericDistriSim, randomizeN0=randomizeN0)
    
    N0<- sum(N0)
    
    if (is.null(distriABM) | all(is.na(distriABM))){
      if (pars$raw){
        rawSim[[n]]<- NA
        names(rawSim)[n]<- paste0("N", N0)
      }
      
      next
    }
    
    
    stats[n, -c(1:2)]<- summary(distriABM) # First columns: idScenario and N0
    
    if (pars$raw){
      rawSim[[n]]<- distriABM
      names(rawSim)[n]<- N0
    }
    if (pars$numericDistriSim){
      distriSim[[n]]<- numericDistriABMSim2numericDistriSim(distriABM)
      names(distriSim)[n]<- N0
    }
    if (pars$Ntf){
      if (pars$numericDistriSim) distriSimTmp<- distriSim[[n]]
      else distriSimTmp<- numericDistriABMSim2numericDistriSim(distriABM)
      
      if (length(distriSimTmp) == 2){ ## Only one timestep or (pars$Ntf & !pars$raw & !pars$numericDistriSim)
        Ntf[[n]]<- distriSimTmp[[2]]
        
      } else {
        tf<- which(sapply(distriSimTmp, inherits, "numericDistri"))
        
        if (length(tf) > 0){
          tf<- tf[length(tf)]
          Ntf[[n]]<- distriSimTmp[[tf]]
        }else{
          Ntf[[n]]<- NA
        }
      }
      
      names(Ntf)[n]<- N0
    }
  }
  
  res<- list(stats=stats)
  if (pars$Ntf) res<- c(res, list(Ntf=Ntf))
  if (pars$numericDistriSim) res<- c(res, list(numericDistriSim=distriSim))
  if (pars$raw) res<- c(res, list(raw=rawSim))
  
  return (res)
}


## result(): Post process result ----
#' @rdname Model
#'
#' @param model 
#' @param type 
#' @details 
#'
#' @return a data frame with the aggregated results and parameters of a simulation.
#' @examples result(res) 
#' @export
setGeneric("result", function(model, type=c("stats", "N0_Pest", "Ntf"), popbio=FALSE, ...) standardGeneric("result"))

setMethod("result", 
          signature(model="Model", type="ANY"),
          function(model, type=c("stats", "N0_Pest", "Ntf"), popbio=FALSE, ...){
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
                res$Pest<- 1 - res$extinct

                rownames(res)<- paste0(res$idScenario, "-N", res$N0)
                res<- res[order(res$idScenario, res$N0),]

                res

              },
              N0_Pest={
                if (nrow(model@sim@N0_Pest) == 0){
                  stop("There are no results yet. Use findN0_Pest(model) to start the simulations. The function return a model with the results on the model@sim@N0_Pest slot.\n")
                }
                res<- merge(S3Part(model), S3Part(model@sim@N0_Pest), by="idScenario")
                
                rownames(res)<- res$idScenario
                res<- res[order(res$idScenario),]
                
                res
                
              },
              Ntf={
                if (inherits(model, "Model.numericDistriABM")){
                  Ntf<- model@sim@Ntf
                  idScenario<- names(Ntf)
                  N0<- as.numeric(sapply(Ntf, names))
                  NtfUL<- unlist(Ntf, recursive=FALSE)
                  
                  NtfUL<- sapply(NtfUL, function(x){
                    quantile(x, na.rm=TRUE, ...)
                  })
                  
                  NtfUL<- t(NtfUL)
                  
                  Ntf<- cbind(data.frame(idScenario, N0, stringsAsFactors=FALSE), NtfUL)
                  
                }else if (inherits(model, "Model.numericDistri")){
                  Ntf<- model@sim@raw
                  Ntf<- unlist(Ntf, recursive=FALSE)
                  
                  Ntf<- sapply(Ntf, function(x){
                    quantile(x, na.rm=TRUE, ...)
                  })
                  
                  Ntf<- t(Ntf)
                  idScenario<- unlist(gsub("\\.[0-9]+$", "", rownames(Ntf)))
                  N0<- as.numeric(unlist(gsub(".+Env[0-9]+\\.", "", rownames(Ntf))))
                  
                  Ntf<- cbind(data.frame(idScenario, N0, stringsAsFactors=FALSE), Ntf)

                }else{
                  
                  if (nrow(model@sim@Ntf) == 0){
                    stop("There are no results yet. Use run(model) to start the simulations. Check that model@sim@params@Ntf is TRUE before running. The function return a model with the results on the @sim@Ntf slot.\n")
                  }
                  
                  Ntf<- apply(model@sim@Ntf[,-c(1,2)], 1, function(x){
                          if (all(is.na(x))){
                            res<- rep(NA_real_, 5) # TODO: add probs parameter for Ntf -> ifelse(missing(probs), rep(NA_real_, 5), rep(NA_real_, length(probs))) ## 5 quantiles otherwise as many as probs length
                          } else {
                            res<- stats::ecdf(x)
                            res<- stats::quantile(res, ...)
                          }
                          res
                        })
                  Ntf<- t(Ntf)
                  colnames(Ntf)<- names(stats::quantile(1, ...))
                  Ntf<- cbind(model@sim@Ntf[, c(1,2)], Ntf, stringsAsFactors=FALSE)
                }
                
                res<- merge(S3Part(model), Ntf, by="idScenario")
                
                rownames(res)<- paste0(res$idScenario, "-N", res$N0)
                
                # Sort rows
                res<- res[order(res$idScenario, res$N0),]
                
                # resl<- reshape2::melt(res[, c("idScenario", "N0", selQuantile)], id.vars=1:2, value.name="Ntf") # id vars: idScenario and N0
                # resl$quantile<- as.numeric(as.character(resl$variable))
                # resl$quantile<- (resl$quantile - 1) / (model@sim@params$replicates - 1) # length(unique(resl$quantile)) == replicates
                # resl<- resl[,-3] # remove "variable" column, a id for replicates
                
                res
              
              }
            )
            
            ## Add colorLH
            idLH<- unique(res$idLH)
            if (any(grepl("-", idLH))){
              lh<- strsplit(idLH, "-")
              lh<- sapply(lh, function(y) y[[1]])
              lh<- data.frame(idLH=idLH, colorLH=lh, stringsAsFactors=FALSE)
              res<- merge(res, lh, by="idLH")
            } else if (length(idLH) < 9){ # default palette() has 8 colors
              res$colorLH<- factor(res$idLH)
            } else {
              res$colorLH<- "black"
            }
            
            if (popbio & requireNamespace("popbio", quietly=TRUE)){
              lh<- unique(res[, c("idLH", "a", "s", "j", "fecundity", "AFR")])
              
              popbio<- apply(lh[,-1], 1, function(x){ # First column is a character and makes x a character vector
                mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=fecundity * j, AFR=AFR))
                return(eigen.analisys2df(mat))
              })
              popbio<- do.call(rbind, popbio)
              popbio<- cbind(idLH=lh$idLH, popbio)
              
              res<- merge(res, popbio, by="idLH")
            }
            
            # sort columns
            res<- cbind(res[, grep("^(idScenario|N0)$", names(res)), drop=FALSE],
                        res[,-c(grep("^(idScenario|N0)$", names(res)))])
            
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
              cat("There are no results yet. Use run(model) to start the simulations for all scenarios and N0 or findN0_Pest(model) to find the N0 giving a fixed probability of extinction.\n")
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
## TODO: use match to keep the index orders from the user input. -> ERRORS in tests
#' @rdname Model
#' @export
`[.Model`<- function(x, ...){
  xSel<- data.frame(x, stringsAsFactors=FALSE)[...]
  scenarioIdSel<- xSel$idScenario
  
  sim<- x@sim
  
  if (nrow(sim) > 0){
    stats<- data.frame(sim)
    stats<- stats[stats$idScenario %in% scenarioIdSel,]
    # stats<- stats[na.omit(match(stats$idScenario, scenarioIdSel)),]
    S3Part(sim)<- stats
  }
  
  if (length(sim@raw) > 0){
    sim@raw<- sim@raw[which(names(sim@raw) %in% scenarioIdSel)]
    # sim@raw<- sim@raw[na.omit(match(names(sim@raw), scenarioIdSel))]
  }
  
  if (nrow(sim@N0_Pest) > 0){
    N0_Pest<- sim@N0_Pest
    sim@N0_Pest<- N0_Pest[N0_Pest$idScenario %in% scenarioIdSel,]
    # sim@N0_Pest<- N0_Pest[na.omit(match(N0_Pest$idScenario, scenarioIdSel)),]
  }
  
  
  if (inherits(sim, "Sim.discretePopSim")){
    if (nrow(sim@Ntf) > 0){
      Ntf<- sim@Ntf
      sim@Ntf<- Ntf[Ntf$idScenario %in% scenarioIdSel,]
      # sim@Ntf<- Ntf[na.omit(match(Ntf$idScenario, scenarioIdSel)),]
    }
    
    if (length(sim@discretePopSim) > 0){
      pop<- sim@discretePopSim
      sim@discretePopSim<- pop[names(pop) %in% scenarioIdSel]
      # sim@discretePopSim<- pop[na.omit(match(names(pop), scenarioIdSel))]
    }
    
  }else if (inherits(sim, "Sim.numericDistri")){
    if (length(sim@Ntf) > 0){
      Ntf<- sim@Ntf
      sim@Ntf<- Ntf[names(Ntf) %in% scenarioIdSel]
      # sim@Ntf<- Ntf[na.omit(match(names(Ntf), scenarioIdSel))]
    }
    
    if (length(sim@numericDistriSim) > 0){
      distriSim<- sim@numericDistriSim
      sim@numericDistriSim<- distriSim[names(distriSim) %in% scenarioIdSel]
      # sim@numericDistriSim<- distriSim[na.omit(match(names(distriSim), scenarioIdSel))]
    }
    
  }

  out<- Model(pars=xSel, sim=sim)
  
  return(out)
}


#' @rdname Model
#' @export
rbind.Model<- function(...){
  models<- list(...)
  
  classModels<- sapply(models, class)
  if (length(unique(classModels)) > 1){
    stop("All model objects must inherit from the same class. There are models inherithing ",
         paste(unique(classModels), collapse=" & "), ".")
  }

  sims<- lapply(models, function(x) x@sim)
  params<- lapply(sims, function(x) x@params)
  params<- unique(params)
  if (length(params) > 1){
    tmpParams<- unique(lapply(params, function(x) x[names(x) != "transitionsFunc"]))
    tmpTransitionsFunc<- lapply(params, function(x) x[names(x) == "transitionsFunc"])
    equalTrans<- sapply(tmpTransitionsFunc, function(x) all.equal(x, tmpTransitionsFunc[[1]]))
    if (length(tmpParams) == 1 & all(equalTrans)){
      params<- params[[1]]
    }else{
      stop("Ensambling models with different simulation parameters is not supported.")
    }
  }else{
    params<- params[[1]]
  }
  
  scenarios<- lapply(models, S3Part)
  idScenarios<- lapply(scenarios, function(x) x$idScenario)
  scenario<- do.call(rbind, scenarios)
  rebuildId<- FALSE
  
  if (any(duplicated(unlist(idScenarios)))){
    nonUniqueIdLH<- FALSE
    nonUniqueIdEnv<- FALSE
    rebuildId<- TRUE
    
    ## duplicated idLH (same idLH & different params)
    colsLH<- c("baseLH", "idLH", "lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR")
    colsLH<- intersect(colsLH, names(scenario))
    colsLHnoId<- colsLH[!grepl("LH", colsLH)] # "baseLH", "idLH"
    
    lhU<- unique(scenario[, colsLHnoId])
    lhUId<- unique(scenario$idLH)
    if (nrow(lhU) != length(lhUId)){
      nonUniqueIdLH<- TRUE
      # lhU$idLH_ori<- lhU$idLH
      lhU<- lhU[order(lhU$lambda, lhU$a, lhU$fecundity), ]
      lhU$idLH<- formatC(seq_len(nrow(lhU)), format="d", flag="0", width=nchar(nrow(lhU)))
    }
    
    
    ## duplicated idEnv
    colsEnv<- c("idEnv", "seasonAmplitude", "seasonMean", "varJ", "varA", "breedFail")
    colsEnv<- intersect(colsEnv, names(scenario))
    colsEnvnoId<- colsEnv[!grepl("idEnv", colsEnv)]
    
    envU<- unique(scenario[, colsEnvnoId])
    envUId<- unique(scenario$idEnv)
    if (nrow(envU) != length(envUId)){
      nonUniqueIdEnv<- TRUE
      # envU$idEnv_ori<- envU$idEnv
      envU<- envU[order(envU$seasonAmplitude, envU$seasonMean, envU$varJ, envU$varA, envU$breedFail), ]
      envU$idEnv<- formatC(seq_len(nrow(envU)), format="d", flag="0", width=nchar(nrow(envU)))
    }
    
    
    ## Rebuild idLH
    if (nonUniqueIdLH){
      scenarios<- lapply(scenarios, function(x){
        res<- merge(x[, -grep("LH", names(x))], lhU, by.x=colsLHnoId, by.y=colsLHnoId)
      })
    }
    
    ## Rebuild idEnv
    if (nonUniqueIdEnv){
      scenarios<- lapply(scenarios, function(x){
        res<- merge(x[, -grep("idEnv", names(x))], envU, by.x=colsEnvnoId, by.y=colsEnvnoId)
      })
    }
    
    
    ## Rebuild idScenario
    if (rebuildId){
      scenarios<- lapply(scenarios, function(x){
        id<- paste0("LH", x$idLH, "_Env", x$idEnv)
        idExtra<- gsub("^LH.+_Env[0-9a-zA-Z]+[_]*", "", x$idScenario)
        if (any(idExtra != "")){
          id<- paste0(id, "_", idExtra)
        }
        
        if (anyDuplicated(id)){
          warning("Some duplicated idScenario values. Be careful and try to fix it before proceed!")
        }
        
        x$idScenario_ori<- x$idScenario
        x$idScenario<- id
        x
      })
    }
    
    colsOrd<- c(names(scenario), setdiff(names(scenarios[[1]]), names(scenario)))
    scenario<- do.call(rbind, scenarios)[, colsOrd]
    
    if (any(dup<- duplicated(scenario$idScenario))){
      dupId<- scenario$idScenario[dup]
      dups<- lapply(scenarios, function(x) x$idScenario_ori[x$idScenario %in% dupId])
      stop("Some duplicated scenarios with the same parameters. Autoremove duplicates not implemented yet. Do it.\n", dupId)
      ## TODO: remove duplicated scenarios keeping the first only and throwing a warning
      warning("Some duplicated scenarios with the same parameters. Remove them.")
      i<- length(scenarios)
      while (length(dupId) > 0){
        m<- models[[i]]
        models[[i]]<- m[!m$idScenario %in% dupId, ]
        scenario<- do.call(rbind, scenarios)[, colsOrd]
        dupId<- scenario$idScenario[duplicated(scenario$idScenario)]
        i<- i - 1
      }
    }
    
    rownames(scenario)<- scenario$idScenario
  }
  
  
  sim<- Sim(params=params, type=gsub("Sim\\.", "", class(sims[[1]])))

  # If nrow == 0, data.frame change from S4 to S3 and some tests fail
  sel<- which(sapply(sims, function(x) nrow(S3Part(x)) > 0))
  if (length(sel) > 0){
    statsL<- lapply(sims[sel], S3Part)
    
    if (rebuildId){
      colsOrd<- names(statsL[[1]])
      
      statsL<- mapply(function(x, y){
                  names(x)<- gsub("idScenario", "idScenario_ori", names(x))
                  res<- merge(x, y[, c("idScenario", "idScenario_ori")], by="idScenario_ori")
                  rownames(res)<- paste0(res$idScenario, "_N", res$N0)
                  res<- res[, colsOrd]
                }, x=statsL, y=scenarios[sel], SIMPLIFY=FALSE)
    }
    
    stats<- do.call(rbind, statsL)
    S3Part(sim)<- stats
  }
  
  sel<- which(sapply(sims, function(x) nrow(x@N0_Pest) > 0))
  if (length(sel) > 0){
    N0_PestL<- lapply(sims[sel], function(x) x@N0_Pest)
    
    if (rebuildId){
      colsOrd<- names(N0_PestL[[1]])
      
      N0_PestL<- mapply(function(x, y){
                    names(x)<- gsub("idScenario", "idScenario_ori", names(x))
                    res<- merge(x, y[, c("idScenario", "idScenario_ori")], by="idScenario_ori")
                    rownames(res)<- res$idScenario
                    res<- res[, colsOrd]
                  }, x=N0_PestL, y=scenarios[sel], SIMPLIFY=FALSE)
    }
    
    sim@N0_Pest<- do.call(rbind, N0_PestL)
  }
  
  
  rawL<- lapply(sims, function(x) x@raw)
  
  if (rebuildId){
    rawL<- mapply(function(x, y){
      if (is.null(x)) return(x)
      
      thesaurus<- merge(y[, c("idScenario", "idScenario_ori")], names(x), by.x="idScenario_ori", by.y=1)
      names(x)<- thesaurus$idScenario[match(names(x), thesaurus$idScenario_ori)]
      x
    }, x=rawL, y=scenarios, SIMPLIFY=FALSE)
  }
  sim@raw<- do.call(c, rawL)

  
  if (inherits(sim, "Sim.discretePopSim")){
    sel<- which(sapply(sims, function(x) nrow(x@Ntf) > 0))
    if (length(sel) > 0){
      NtfL<- lapply(sims[sel], function(x) x@Ntf)
      
      if (rebuildId){
        colsOrd<- names(NtfL[[1]])
        
        NtfL<- mapply(function(x, y){
          names(x)<- gsub("idScenario", "idScenario_ori", names(x))
          res<- merge(x, y[, c("idScenario", "idScenario_ori")], by="idScenario_ori")
          rownames(res)<- paste0(res$idScenario, "_N", res$N0)
          res<- res[, colsOrd]
        }, x=NtfL, y=scenarios[sel], SIMPLIFY=FALSE)
      }
      
      sim@Ntf<- do.call(rbind, NtfL)
    }
    
    
    discretePopSimL<- lapply(sims, function(x) x@discretePopSim)
    
    if (rebuildId){
      discretePopSimL<- mapply(function(x, y){
        if (is.null(x) | length(x) == 0) return(x)
        
        thesaurus<- merge(y[, c("idScenario", "idScenario_ori")], names(x), by.x="idScenario_ori", by.y=1)
        names(x)<- thesaurus$idScenario[match(names(x), thesaurus$idScenario_ori)]
        x
      }, x=discretePopSimL, y=scenarios, SIMPLIFY=FALSE)
    }
    
    sim@discretePopSim<- do.call(c, discretePopSimL)
  }
  
  if (inherits(sim, "Sim.numericDistri")){
    sel<- which(sapply(sims, function(x) length(x@Ntf) > 0))
    if (length(sel) > 0){
      NtfL<- lapply(sims[sel], function(x) x@Ntf)
      
      if (rebuildId){
        NtfL<- mapply(function(x, y){
          if (is.null(x) | length(x) == 0) return(x)
          
          thesaurus<- merge(y[, c("idScenario", "idScenario_ori")], names(x), by.x="idScenario_ori", by.y=1)
          names(x)<- thesaurus$idScenario[match(names(x), thesaurus$idScenario_ori)]
          x
        }, x=NtfL, y=scenarios, SIMPLIFY=FALSE)
      }
      
      sim@Ntf<- do.call(c, NtfL)
    }
    
    
    distriSimL<- lapply(sims, function(x) x@numericDistriSim)
    
    if (rebuildId){
      distriSimL<- mapply(function(x, y){
        if (is.null(x) | length(x) == 0) return(x)
        
        thesaurus<- merge(y[, c("idScenario", "idScenario_ori")], names(x), by.x="idScenario_ori", by.y=1)
        names(x)<- thesaurus$idScenario[match(names(x), thesaurus$idScenario_ori)]
        x
      }, x=distriSimL, y=scenarios, SIMPLIFY=FALSE)
    }
    
    sim@numericDistriSim<- do.call(c, distriSimL)
  }
  
  
  return(Model(pars=scenario, sim=sim))
}


## TODO: numericDistri plots ----
#' Plot Model
#'
#' @rdname Model
#' @param x 
#' @param resultType 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom graphics plot
plot.Model<- function(x, resultType=c("Pest_N0", "G", "N0_Pest", "Ntf"), ...){
  
  if (missing(resultType)) noType<- TRUE else noType<- FALSE
  
  resultType<- match.arg(resultType)
  
  stats<- nrow(x@sim) > 0
  N0_Pest<- nrow(x@sim@N0_Pest) > 0
  if (inherits(x, "Model.numericDistri")){
    Ntf<- length(x@sim@raw) > 0
  }else{
    Ntf<- nrow(x@sim@Ntf) > 0
  }
  
  # if no type is specified select a existing one with precedence stats > N0_Pest > Ntf
  if (noType){
    if (stats) resultType<- "Pest_N0"
    else if (N0_Pest) resultType<- "N0_Pest"
    else if (Ntf) resultType<- "Ntf"
    else {
      message("No results found. Plotting the parameter space.\n\trun(model) or findN0_Pest(model) to simulate.")
      x<- S3Part(x)
      if ("baseLH" %in% names(x)){
        x$colorLH<- factor(x$baseLH)
      } else {
        x$colorLH<- "black"
      }
      
      cols<- intersect(names(x), c("lambda", "fecundity", "broods", "b", "a", "s", "j", "AFR",
                "interval", "seasonAmplitude", "seasonMean", "varJ", "varA", "breedFail", "jind", "jbr", "colorLH"))
      x<- unique(x[, cols])
      out<- graphics::plot(x[, sapply(x, is.numeric)], col=x$colorLH, ...) # All selected columns except colorLH
      # graphics::legend("topright", legend=levels(res$colorLH), bty = "y", pch = 19, col=res$colorLH)
      
      return(invisible(out))
    }
  }
  
  out<- NA
  
  if (stats & resultType == "Pest_N0"){
    res<- result(x, type="stats")

    out<- plotPest_N0(res, ...)
  }
  
  if (stats & resultType == "G"){
    res<- result(x, type="stats")

    out<- plotG(res, ...)
  }
  
  if (N0_Pest & resultType == "N0_Pest"){  
    res<- result(x, type="N0_Pest")

    out<- plotN0_Pest(res, ...)
  }

  if (Ntf & resultType == "Ntf"){
    res<- result(x, type="Ntf")

    out<- plotNtf(res, ...)
  }
  
  # out<- out + ggplot2::scale_color_gradient2(low="red", mid="yellow", high="green") +
  #   ggplot2::scale_fill_gradient2(low="red", mid="yellow", high="green")
  
  # out<- out + ggplot2::scale_color_gradientn(colors=grDevices::terrain.colors(nrow(lh))) +
  #   ggplot2::scale_fill_gradientn(colors=grDevices::terrain.colors(nrow(lh)))
  
  return(out)
}


plotPest_N0<- function(x, ...){
  x$Pest<- 1 - x$extinct

  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, y=Pest, group=idScenario, color=colorLH)) + 
    ggplot2::geom_line(mapping=ggplot2::aes(size=lambda), alpha=0.5) + ggplot2::geom_point() +
    ggplot2::facet_grid(seasonAmplitude + breedFail ~ varA + varJ, labeller=ggplot2::label_both) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2))
}

plotG<- function(x, ...){
  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, y=GL, group=idScenario, color=colorLH)) + 
    ggplot2::geom_hline(yintercept=1) + ggplot2::geom_line(mapping=ggplot2::aes(size=lambda), alpha=0.5) + # ggplot2::geom_point() +
    ggplot2::geom_line(mapping=ggplot2::aes(y=meanL, size=lambda), linetype=2) +
    ggplot2::facet_grid(seasonAmplitude + breedFail ~ varA + varJ, labeller=ggplot2::label_both) +
    ggplot2::labs(x=expression(N[0]), y=expression(lambda * " geometric mean and mean (dashed) for all transitions")) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2))
}

plotN0_Pest<- function(x, ...){
  ggplot2::ggplot(data=x, ggplot2::aes(x=lambda, y=N0interpoled, group=colorLH, color=colorLH)) + ggplot2::scale_y_log10() +
    ggplot2::geom_point() + ggplot2::geom_line(size=0.5) + ggplot2::facet_grid(seasonAmplitude + breedFail ~ varA + varJ, labeller=ggplot2::label_both)  
}

plotNtf<- function(x, ...){
  abline<- 1:max(x$N0)
  abline<- data.frame(x=abline, y=abline)
  
  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, color=colorLH, fill=colorLH, group=idScenario)) + ggplot2::scale_y_log10() + 
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=`25%`, ymax=`75%`), alpha=.4, linetype=0) + ggplot2::geom_line(mapping=ggplot2::aes(y=`50%`, size=lambda)) + 
    ggplot2::facet_grid(seasonAmplitude + breedFail ~ varA + varJ, labeller=ggplot2::label_both) + ggplot2::labs(x=expression(N[0]), y=expression(N[tf] * " (mean and 50%)")) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2)) + ggplot2::geom_line(mapping=ggplot2::aes(x, y), data=abline, inherit.aes=FALSE, show.legend=FALSE, lty=2)
}


#' Plot a histogram with the final population size of each replicate
#'
#' @rdname Model
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom graphics hist
#' @examples
hist.Model<- function(x, ...){
  if (nrow(x@sim) == 0){
    stop("No results found.\n\trun(model) to simulate.")
  }
  
  Ntf<- x@sim@Ntf
  xd<- data.frame(x, stringsAsFactors=FALSE)
  
  Ntf<- merge(Ntf, xd[, c("idScenario", "seasonAmplitude", "varJ", "varA", "breedFail")], by="idScenario")
  N<- reshape2::melt(Ntf, id.vars=c("idScenario", "N0", "seasonAmplitude", "varJ", "varA", "breedFail"), value.name="Ntf")
  
  ggplot2::ggplot(N, ggplot2::aes(x=Ntf, group=idScenario, color=idScenario)) + 
    # ggplot2::geom_histogram() + 
    ggplot2::geom_freqpoly(show.legend=FALSE) + # ggplot2::scale_x_log10() +
    ggplot2::facet_grid(seasonAmplitude + breedFail ~ N0 + varA + varJ, labeller=ggplot2::label_both)
}
