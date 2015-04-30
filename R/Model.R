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
## Parametres ----
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


