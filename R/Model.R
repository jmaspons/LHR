## Temporal variability class
setClass("Model", slots=list(sim="Sim", params="list"), contains="data.frame")

## Constructor ----
setGeneric("Model", function(lh=LH(), env=Env(), sim=Sim()) standardGeneric("Model"))

setMethod("Model",
          signature(lh="missing", env="missing", sim="missing"),
          function(lh=LH(), env=Env(), sim=Sim()){
            model<- Model(lh, env, sim)
            return (model)
          }
)


setMethod("Model",
          signature(lh="LH", env="Env", sim="Sim"),
          function(lh, env, sim){
            lhEnv<- combineLH_Env(lh=lh, env=env)
            scenario<- lhEnv$scenario
            parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv, breedFail=lhEnv$breedFail)
            
            model<- new("Model", 
                        scenario,
#                         lh=lh, Data in scenario and easily to recreate the obj with LH(unique(scenario))
#                         env=env,
                        sim=sim,
                        params=parameters)
            return (model)
          }
)

setMethod("Model",
          signature(lh="LH", env="Env", sim="missing"),
          function(lh, env){
            Model(lh=lh, env=env, sim=Sim())
          }
)


## Combine LH and Environment ----
setGeneric("combineLH_Env", function(lh=LH(), env=Env(), resolution=12, nSteps, interval=2, criterion=c("maxFirst", "maxMean")[1]) standardGeneric("combineLH_Env"))

# return a sinusoidal pattern
setMethod("combineLH_Env", 
          signature(lh="missing", env="missing", resolution="missing", nSteps="missing", interval="missing", criterion="missing"),
          function(lh=LH(), env=Env()){
            combineLH_Env(lh, env)
          }
)

setMethod("combineLH_Env", 
          signature(lh="LH", env="Env", resolution="ANY", nSteps="ANY", interval="ANY", criterion="ANY"),
          function(lh, env, resolution=12, nSteps, interval=1, criterion="maxFirst"){
            scenario<- merge(S3Part(lh), S3Part(env))
            
            ## Brood failure (e.g. nest predation)
            # breed failure is a proportion of juvenile mortality correlated at the brood level j = jEnv * breedFailLH
            #  j is survival, and breed failure a mortality
            scenario$jM<- 1 - (1 - scenario$j) * (1 - scenario$breedFail)  # juvenile mortality is divided in brood mortality and individual mortality
            scenario$breedFailM<- (1 - scenario$j) * scenario$breedFail

#             sel<- which(scenario$breedFail == 1)
#             scenario$jM[sel]<- 1 # corrects /0.  When breedFail == 1 all juvenile mortality is from breeding fail and there is no individual mortality
#             scenario$breedFailM<- scenario$j * scenario$breedFail
            
            ## Seasonality
            seasonBroods<- data.frame(unique(scenario[,c("mean", "seasonAmplitude", "broods")]), var=0, breedFail=0) # var necessary for Env(seasonBroods)
            seasonVar<- seasonOptimCal(env=Env(seasonBroods), resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
            
            ##Stochasticity
            # Calculated later when running simulations
#             betaPars<- unique(scenario[,c("mean", "var", "a", "s", "j", "jEnv")])
#             fbeta(mean=betaPars$mean * betaPars, var=scevar)
            
            
            return (list(scenario=scenario, seasonBroodEnv=list(parsBroodEnv=seasonBroods, broodEnv=seasonVar))) #, breedFail=breedFail))
          }
)

# returns the nSteps values separed by interval units of a sinusoidal pattern optimizing 
# according to different criterion (maxMean, maxFirst)
# seasonEvents<- function(env)

## Simulate ----
setGeneric("run", function(model, ...) standardGeneric("run"))

# Create a Sim object with the results
setMethod("run", 
          signature(model="Model"),
          function(model, ...){
            simRes<- switch(class(model@sim),
                               Sim.discretePopSim=run.discretePopSim(model),
                               Sim.numericDistri=run.numericDistri(model),
                               Sim.ssa=run.ssa(model, ...))
# print(str(simRes))
            modelRes<- new("Model", 
                        S3Part(model),
                        sim=simRes,
                        params=model@params)
            return(modelRes)
          }
          
)

run.discretePopSim<- function(model){
# print("Inside run.discretePop()")
  scenario<- S3Part(model)
  pars<- model@sim@params

  rawSim<- list()
  res<- matrix(nrow=nrow(scenario) * length(pars$N0), ncol=12, 
              dimnames=list(scenario=paste0(rep(rownames(scenario), each=length(pars$N0)), "_N", rep(pars$N0, times=nrow(scenario))),
                            stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  if (pars$Ntf){
    Ntf<- matrix(nrow=nrow(scenario) * length(pars$N0), ncol=2 + pars$replicates, 
                 dimnames=list(scenario=paste0(rep(rownames(scenario), each=length(pars$N0)), "_N", rep(pars$N0, times=nrow(scenario))), 
                               Ntf=c("scenario", "N0", 1:pars$replicates)))
    Ntf<- as.data.frame(Ntf)
  }
  
  k<- 1
  for (i in 1:nrow(scenario)){ ## TODO: parallelise
    ## TODO: pass model@params (LH*Env) to the function, add as a scenario instead as in @params
    rawSim[[i]]<- list()
    ## Seasonality
    seasonVar<- seasonOptimCal(env=Env(scenario[i,]))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
    if (any(seasonVar !=1)){
      jMseason<- scenario$jM[i] * seasonVar
      breedFailMseason<- scenario$breedFailM[i] * seasonVar
    }else{
      jMseason<- scenario$jM[i]
      breedFailMseason<- scenario$breedFailM[i]
    }
    
    for (n in 1:length(pars$N0)){
      N0<- pars$N0[n]
      pop<- with(scenario[i,], discretePopSim(broods=broods, b=b, j=jMseason, a=a, breedFail=breedFailMseason,
                                              varJ=ifelse(pars$envVar$j, var, 0), varBreedFail=ifelse(pars$envVar$breedFail, var, 0),
                                              pars$sexRatio, pars$matingSystem, N0, pars$replicates, pars$tf))
      ## TODO: pop<- list(popF, popM)
      res[k,]<- c(i, N0, as.numeric(summary(pop)))

      if (pars$raw){
        rawSim[[i]][[n]]<- pop
        names(rawSim[[i]])[n]<- N0
      }
      if (pars$Ntf){
        pop<- pop[,ncol(pop)]
        pop[is.na(pop)]<- 0
        Ntf[k,]<- c(rownames(scenario)[i], N0, sort(pop))
      }

      k<- k + 1
    }
  }
  names(rawSim)<- row.names(scenario)

  res<- as.data.frame(res)
  simRes<- new("Sim.discretePopSim", res, params=pars)
  if (pars$raw){
    simRes@raw<- rawSim
  }
  if (pars$Ntf){
    Ntf[,-1]<- apply(Ntf[,2:ncol(Ntf)], 2, as.numeric)
    simRes@Ntf<- Ntf
  }
  
  return(simRes)
}

run.numericDistri<- function(model){
  # print("Inside run.numericDistri()")
  scenario<- S3Part(model)
  pars<- model@sim@params
  
  rawSim<- list()
  res<- matrix(nrow=nrow(scenario) * length(pars$N0), ncol=12, 
              dimnames=list(scenario=paste0(rep(rownames(scenario), each=length(pars$N0)), "_N", rep(pars$N0, times=nrow(scenario))),
                            stats=c("scenario", "N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  k<- 1
  for (i in 1:nrow(scenario)){ ## TODO: parallelise
    ## TODO: pass model@params (LH*Env) to the function, add as a scenario instead as in @params
    rawSim[[i]]<- list()
    ## Seasonality
    seasonVar<- seasonOptimCal(env=Env(scenario[i,]))[[1]] #, resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)
    if (any(seasonVar !=1)){
      jMseason<- scenario$jM[i] * seasonVar
      breedFailMseason<- scenario$breedFailM[i] * seasonVar
    }else{
      jMseason<- scenario$jM[i]
      breedFailMseason<- scenario$breedFailM[i]
    }
print(scenario[i,])    
    for (n in 1:length(pars$N0)){
      N0<- pars$N0[n]
      distri<- with(scenario[i,], t1distri(broods=broods, b=b, j=jMseason, a=a, breedFail=breedFailMseason,
                                              varJ=ifelse(pars$envVar$j, var, 0), varBreedFail=ifelse(pars$envVar$breedFail, var, 0),
                                              pars$sexRatio, pars$matingSystem, N0))
      
      if (pars$raw){
        rawSim[[i]][[n]]<- distri
        names(rawSim[[i]])[n]<- N0
      }
      ## TODO: pop<- list(popF, popM)
      res[k,]<- c(i, N0, as.numeric(summary(distri)))
      k<- k + 1
    }
  }
  names(rawSim)<- row.names(params)
  
  res<- as.data.frame(res)
  if (pars$raw){
    simRes<- new("Sim.numericDistri", res, params=pars, raw=raw)
  }else{
    simRes<- new("Sim.numericDistri", res, params=pars)
  }
  
  return(simRes)
}

##
setGeneric("result", function(model, type="stats") standardGeneric("result"))
# Create a data frame with the aggregated results and parameters
setMethod("result", 
          signature(model="Model", type="ANY"),
          function(model, type="stats"){
            if (nrow(model@sim) == 0){
              stop("There are no results yet. Use run(model) to start the simulations. The function return a model with the results on the sim slot.\n")
            }
#               if (all(is.na(as.numeric(res$scenario)))){ # for categoric scenarios in ssaLH-behavior
#                 scenario<- strsplit(res$scenario, "_")
#                 scenario<- do.call("rbind", scenario)
#                 if (inherits(model@sim, "Sim.ssa")){
#                   colnames(scenario)<- c("LH", "environment", "behavior")
#                 }
#                 res<- data.frame(scenario, res)
#               }
            res<- switch(type,
              stats={
                res<- merge(S3Part(model), S3Part(model@sim), by.x="row.names", by.y="scenario")
                names(res)[1]<- "scenario"
                res<- res[order(res$scenario, res$N0),]
                res},
              Ntf={
                res<- model@sim@Ntf
                res<- melt(res, id.vars=1:2, value.name="Ntf") # id vars: scenario and N0
                res$quantile<- as.numeric(res$variable)
                res$quantile<- (res$quantile - 1) / (model@sim@params$replicates - 1) # length(unique(res$quantile)) == replicates
                res<- res[,-3]
              }
            )

          return(res)
        }
)


## Generic methods ----
setMethod("print", signature(x="Model"),
          function(x, ...){
            print(S3Part(x), ...) # S3Part(x)[,1:12]
          }
)

setMethod("show", signature(object="Model"),
          function(object){
            cat("Object of class \"Model\" with", nrow(object), "scenarios\n")
            if (nrow(object@sim) == 0){
              print(S3Part(object)) # S3Part(x)[,1:12]
              cat("There are no results yet. Use run(model) to start the simulations.\n")
            }else{
              res<- merge(S3Part(object), S3Part(object@sim), by.x="row.names", by.y="scenario")
              names(res)[1]<- "scenario"
              res<- res[order(res$scenario, res$N0),]
              print(res)
              cat("Use result(model) to get a data.frame with the parameters and the results.\n")
            }
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.Model`<- function(x, ...){
  Model(S3Part(x)[...])
}


