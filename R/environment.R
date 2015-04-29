## Temporal variability
setClass(Environment, slots=list(parameters="data.frame", parametersBeta="data.frame", seasonalRange="numeric"))

Environment<- function(mean=1, var=0, seasonAmplitude, seasonRange){
  if (!missing(seasonAmplitude) & !missing(seasonalRange))
    stop("Incorrect parameters for seasonality. Specify seasonRange or seasonAmplitude but not both.")
  
  if (!missing(seasonAmplitude)){
    seasonRange<- seasonRange(mean, amplitude)
  }
  if (!missing(seasonRange)){
    seasonAmplitude<- seasonAmplitude(seasonRange)
    mean<- seasonAmplitude$mean
    seasonAmplitude<- seasonAmplitude$seasonAmplitude
  }
  
  params<- data.frame(mean=mean, var=var, seasonAmplitude=seasonAmplitude)
  paramsBeta<- ifelse(var == 0, NA, fbeta(mean=mean, var=var))
  
  env<- new("Environment", 
            parameters=params,
            parametersBeta=paramsBeta,
            seasonalRange=as.numeric(seasonRange))
  return (env)
}

setGeneric("seasonalPattern")

setMethod("seasonalPattern",
          c(env="Environment", nSteps="integer", interval="numeric", resolution=12, criterion=c("maxFirst", "maxMean")),
          function(env, nSteps, interval, resolution, criterion){
#             seasonality(1, env@d)
          })

# Seasonality

# Stochasticity

# Temporal correlation
#
