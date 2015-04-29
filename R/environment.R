## Temporal variability class
setClass("Environment", slots=list(parameters="data.frame", parametersBeta="data.frame", seasonRange="numeric"))

## Constructor ----
setGeneric("Environment", function(mean=1, var=0, seasonAmplitude=0, seasonRange) standardGeneric("Environment"))

setMethod("Environment",
          signature(mean="numeric", var="ANY", seasonAmplitude="ANY", seasonRange="missing"),
          function(mean=1, var=0, seasonAmplitude=0){
            params<- data.frame(mean=mean, var=var, seasonAmplitude=seasonAmplitude)
            if(var == 0){
              paramsBeta<- data.frame(shape1=NA, shape2=NA)
            }else{
              paramsBeta<- fbeta(mean=mean, var=var)
            }
            
            env<- new("Environment", 
                      parameters=params,
                      parametersBeta=paramsBeta,
                      seasonRange=seasonRange(mean, seasonAmplitude))
            return (env)
          }
)

setMethod("Environment",
          signature(mean="missing", var="ANY", seasonAmplitude="missing", seasonRange="numeric"),
          function(mean, var=0, seasonRange){
            seasonAmplitude<- seasonPars(seasonRange)
            mean<- seasonAmplitude$mean
            seasonAmplitude<- seasonAmplitude$seasonAmplitude
            
            env<- Environment(mean=mean, var, seasonAmplitude=seasonAmplitude)
            return (env)
          }
)  

## Seasonal pattern ----
setGeneric("seasonalPattern", function(env, resolution=12, cicles=1, nSteps, interval, criterion=c("maxFirst", "maxMean")) standardGeneric("seasonalPattern"))

# return a sinusoidal pattern
setMethod("seasonalPattern", 
          signature(env="Environment", resolution="ANY", cicles="ANY", nSteps="missing", interval="missing", criterion="missing"),
          function(env, resolution=12, cicles=1){ 
            timeSteps<- seq(0, 2*pi*cicles - 2*pi/resolution, by=2*pi/resolution)
            pat<- sin(timeSteps) * env@parameters$seasonAmplitude/2 + env@parameters$mean
            return (pat)
          }
)

# returns the nSteps values separed by interval units of a sinusoidal pattern optimizing 
# according to different criterion (maxMean, maxFirst)
setMethod("seasonalPattern",
          signature(env="Environment", resolution="ANY", cicles="ANY", 
                    nSteps="numeric", interval="numeric", criterion="ANY"),
          function(env, resolution=12, cicles=1, nSteps, interval, criterion=c("maxFirst", "maxMean")){
            pattern<- seasonalPattern(env, resolution, cicles)
            n<- length(pattern)
            if (n < nSteps * interval) stop("Resolution is too low for the number of steps and the interval parameters")
            
            if (criterion[1] == "maxFirst"){
              firstDate<- which.max(pattern)
              dates<- seq(firstDate, firstDate + nSteps * interval, length=nSteps)
              if (firstDate + nSteps * interval > n){
                nextPeriod<- which(firstDate + 1:nSteps * interval > n)
                dates[nextPeriod]<- dates[nextPeriod] - n
              }
              
            }else if (criterion[1] == "maxMean"){
              seasonLength<- nSteps * interval
              dates<- seq(1, seasonLength, by=interval)
              ##TODO Optimization
              #     ?optimize(function(i, pattern, dates) sum(pattern[dates + round(i)]), interval=0:n, maximum=TRUE tol=1)
              pattern<- rep(pattern,2)
              objective<- -Inf
              for (i in 0:n){
                newObjective<- sum(pattern[dates + i])
                if (objective < newObjective){
                  objective<- newObjective
                  bestFirstDate<- i
                }
              }
              dates<- dates + bestFirstDate
            }
            return(pattern[dates])
          }
)

# Seasonality auxiliary functions ----
# resolution: time steps in one cicle (e.g. 12 for months, 365 for days)
seasonPars<- function(seasonRange){
  mean<- (min(seasonRange) + max(seasonRange))/2
  seasonAmplitude<- max(seasonRange) - min(seasonRange)
  
  return (data.frame(mean, seasonAmplitude))
}

seasonRange<- function(mean, seasonAmplitude){
  return (c(min=mean - seasonAmplitude/2, max=mean + seasonAmplitude/2))
}

## Stochasticity
# beta distributded value

## Temporal correlation
# TODO
