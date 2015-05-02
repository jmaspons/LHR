# TODO: Add nest mortality. Additive or a part of j??
  #   LH$j<- LH$j / (1 - nestFail)
  #   to-report annual2monthMortality [p]      
  #   annualSurv = monthlySurv ^ 12    =>   monthlySurv = anualSurv ^ 1 / 12
  #   monthSurv = (1 - annualMortality) ^ (1 / 12)
  #   negative binomial: annualSurv = pnbinom(q=0, size=12, prob=monthlySurv) => 0 fails en 12 trails
  #   report 1 - (1 - p) ^ (1 / 12)

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
            parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv)
            
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
            
#             if (missing(resolution)) resolution<- 12 # months
#             if (missing(nSteps)) nSteps<- scenario$broods
#             if (missing(nSteps)) interval<- 1 # month
#             if (missing(criterion)) criterion<- "maxFirst"
            
            ## Seasonality
            seasonBroods<- data.frame(unique(scenario[,c("mean", "seasonAmplitude", "broods")]), var=0)
            seasonVar<- seasonOptimCal(env=Env(seasonBroods), resolution=resolution, nSteps=seasonBroods$broods, interval=interval, criterion=criterion)

            return (list(scenario=scenario, seasonBroodEnv=list(parsBroodEnv=seasonBroods, broodEnv=seasonVar)))
          }
)

# returns the nSteps values separed by interval units of a sinusoidal pattern optimizing 
# according to different criterion (maxMean, maxFirst)
# seasonEvents<- function(env)

## Simulate ----
setGeneric("run", function(model) standardGeneric("run"))

# Create a Sim object with the results
setMethod("run", 
          signature(model="Model"),
          function(model){
            simRes<- switch(class(model@sim),
                               Sim.discretePopSim=run.discretePopSim(model),
                               Sim.numericDistri=run.numericDistri(model),
                               Sim.ssa=run.ssa(model))
# print(str(simRes))
            modelRes<- new("Model", 
                        S3Part(model),
                        sim=simRes,
                        params=model@params)
print(str(modelRes))
            return(modelRes)
          }
          
)

run.discretePopSim<- function(model){
print("Inside run.discretePop()")
  pars<- model@sim@params
  scenario<- S3Part(model)
  rawSim<- list()
  res<- array(dim=c(nrow(scenario) * length(pars$N0), 11), 
              dimnames=list(LH_N0=NULL, stats=c("N0", "increase", "decrease", "stable", "extinct", "GR", "meanR", "varR", "GL", "meanL", "varL"))) # 11 = ncol(summary(pop)) + 1 (N0)
  k<- 1
  for (n in 1:length(pars$N0)){
    N0<- pars$N0[n]
    rawSim[[n]]<- list()
    for (i in 1:nrow(scenario)){ ## TODO: parallelise
      pop<- switch(pars$structure,
              fit=with(scenario[i,], mFit.t(fecundity=fecundity, j=j, a=a, N0=N0, replicates=pars$replicates, tf=pars$tf)),
              survBV=with(scenario[i,], mSurvBV.t(broods=broods, b=b, j=j, a=a, nestFail, N0=N0, replicates=pars$replicates, tf=pars$tf)),
              fitSex=with(scenario[i,], mFitSex.t(fecundity=fecundity, j=j, a=a, sexRatio, matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf)),
              survBVSex=with(scenario[i,], mSurvBVSex.t(broods=broods, b=b, j=j, a=a, nestFail, sexRatio, matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf)))
      if (pars$raw){
        rawSim[[n]][[i]]<- pop
      }
      res[k,]<- c(N0, as.numeric(summary(pop)))
      k<- k + 1
    }
  }
  names(rawSim)<- pars$N0

  res<- as.data.frame(res)
  simRes<- new("Sim.discretePopSim", res, params=pars)
  
  return(simRes)
}

run.numericDistr<- function(){}
run.ssa<- function(){}

## Generic methods ----
setMethod("print", signature(x="Model"),
          function(x, ...){
            print(S3Part(x)[,1:12], ...)
          }
)

setMethod("show", signature(object="Model"),
          function(object){
            cat("Object of class \"Model\" with", nrow(object), "scenarios\n")
            print(S3Part(object)[,1:11])
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns
`[.Model`<- function(x, ...){
  Model(S3Part(x)[...])
}


