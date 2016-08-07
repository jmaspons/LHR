## Temporal variability class

#' Environment class
#'
#' @slot seasonRange data.frame. 
#'
#' @examples Env()
#' @export
setClass("Env", slots=list(seasonRange="data.frame"), contains="data.frame")

## Constructor ----

#' @describeIn Env
#'
#' @param mean 
#' @param var 
#' @param seasonAmplitude 
#' @param seasonRange 
#' @param breedFail 
#'
#' @return a Env class object
#' @export
setGeneric("Env", function(mean=1, var=0, seasonAmplitude=0, seasonRange, breedFail=0) standardGeneric("Env"))

setMethod("Env",
          signature(mean="ANY", var="ANY", seasonAmplitude="ANY", seasonRange="missing", breedFail="ANY"),
          function(mean=1, var=0, seasonAmplitude=0, breedFail=0){
            params<- data.frame(mean=mean, var=var, seasonAmplitude=seasonAmplitude, breedFail=breedFail)#, fbeta(mean=mean, var=var))
#             params[params$var == 0, c("shape1", "shape2")]<- NA
#             params<- params[!is.nan(params$shape1),] ## Remove combinations of mean and var out of the domain of the Beta distribution
            env<- new("Env", 
                      params,
                      seasonRange=seasonRange(mean, seasonAmplitude))
            return (env)
          }
)

setMethod("Env",
          signature(mean="missing", var="missing", seasonAmplitude="missing", seasonRange="missing", breedFail="missing"),
          function(){
            season<- data.frame(mean=c(1,.5), seasonAmplitude=c(0,1))
            comb<- expand.grid(list(var=c(0,0.1), breedFail=c(0, 0.5, 1)))
            comb<- merge(season, comb)
            env<- with(comb, Env(mean=mean, var=var, seasonAmplitude=seasonAmplitude, breedFail=breedFail))
            return (env)
          }
)

setMethod("Env",
          signature(mean="missing", var="ANY", seasonAmplitude="missing", seasonRange="matrix", breedFail="ANY"),
          function(var=0, seasonRange, breedFail=0){
            mat<- matrix(seasonRange, ncol=2)
            seasonPars<- seasonPars(mat)
            seasonMean<- seasonPars$seasonMean
            seasonAmplitude<- seasonPars$seasonAmplitude
            
            env<- Env(mean=seasonMean, var=var, seasonAmplitude=seasonAmplitude, breedFail=breedFail)
            return (env)
          }
)  

setMethod("Env",
          signature(mean="data.frame", var="missing", seasonAmplitude="missing", seasonRange="missing", breedFail="missing"),
          function(mean){
            if (inherits(mean, "Model")) mean<- S3Part(mean)
            x<- unique(mean[,c("mean", "var", "seasonAmplitude", "breedFail")])
            env<- Env(mean=x$mean, var=x$var, seasonAmplitude=x$seasonAmplitude, breedFail=x$breedFail)
            return (env)
          }
)  

## Seasonal pattern ----
#' @describeIn Env 
#'
#' @param env 
#' @param resolution 
#' @param cicles
#' 
#' @return A numeric vector with a sinusoidal pattern
#'
#' @export
setGeneric("seasonalPattern", function(env, resolution=12, cicles=1) standardGeneric("seasonalPattern"))


setMethod("seasonalPattern", 
          signature(env="Env", resolution="ANY", cicles="ANY"),
          function(env, resolution=12, cicles=1){ 
            timeSteps<- seq(0, 2*pi*cicles - 2*pi/resolution, by=2*pi/resolution)
            pat<- matrix(ncol=length(timeSteps), nrow=nrow(env), dimnames=list(NULL, t=1:length(timeSteps)))
            for (i in 1:nrow(env)){
              pat[i,]<- sin(timeSteps) * env$seasonAmplitude[i]/2 + env$mean[i]
            }
            return (pat)
          }
)

# Align a number of events separed by an interval to maximize a seasonal pattern of Env
#' @describeIn Env
#'
#' @param env 
#' @param resolution 
#' @param nSteps 
#' @param interval 
#' @param criterion 
#'
#' @export
setGeneric("seasonOptimCal", function(env, resolution=12, nSteps=2, interval=1, criterion="maxFirst") standardGeneric("seasonOptimCal"))

#' @export
setMethod("seasonOptimCal", 
          signature(env="Env", resolution="missing", nSteps="missing", interval="missing", criterion="missing"),
          function(env){
            seasonOptimCal(env, resolution=12, nSteps=2, interval=1, criterion="maxFirst")
          }
)

#' @export  
setMethod("seasonOptimCal", 
          signature(env="Env", resolution="numeric", nSteps="numeric", interval="numeric", criterion="character"),
          function(env, resolution, nSteps, interval, criterion){
            pattern<- seasonalPattern(env=env, resolution=resolution, cicles=1)
            n<- ncol(pattern) # resolution
            # Extend patterns 2 cicles to displace the event dates and find the best mean
            pattern<- matrix(rep(pattern,2), ncol=n * 2)
            
            # event calendar. SeasonAmplitude ensures that nrow(e) = nrow(env)
            e<- data.frame(nSteps, interval, seasonLength=nSteps * interval, season=env$seasonAmplitude, stringsAsFactors=FALSE)
            ## Maximize pattern on events
            if (criterion[1] == "maxFirst"){
              firstDate<- apply(pattern, 1, which.max)
              dates<- list()
              for (i in 1:nrow(e)){
                if (n < e$nSteps[i] * e$interval[i]){
                  warning("Time steps don't fit in one cicle because number of steps (",
                          e$nSteps[i], ") and the interval (", e$interval[i], ") spans more than the resolution (", resolution, ") of one cicle.")
                  dates[[i]]<- NA
                }else{
                  dates[[i]]<- seq(firstDate[i], firstDate[i] + e$nSteps[i] * e$interval[i], length=e$nSteps[i])
                }
              }
            }else if (criterion[1] == "maxMean"){
              dates<- list()
              for (i in 1:nrow(e)){
                dates[[i]]<- seq(1, e$seasonLength[i], by=e$interval[i])
                
                ##TODO Optimization
                #     ?optimize(function(i, pattern, dates) sum(pattern[dates + round(i)]), interval=0:n, maximum=TRUE tol=1)
                
                objective<- -Inf
                ## Move dates and find the best mean
                for (j in 0:n){
                  newObjective<- sum(pattern[i,dates[[i]] + j])
                  if (objective < newObjective){
                    objective<- newObjective
                    bestFirstDate<- j
                  }
                }
                dates[[i]]<- dates[[i]] + bestFirstDate
              }
            }
            
            eventVal<- list()
            for (i in seq_along(dates)){
              eventVal[[i]]<- pattern[i, dates[[i]]]
            }
            
            return(eventVal)
          }
)

# setMethod("seasonOptimCal", 
#           signature(env="data.frame", resolution="missing", nSteps="missing", interval="missing", criterion="missing"),
#           function(env){
#             seasonOptimCal(Env(env), resolution=env$resolution, nSteps=env$nSteps, interval=env$interval, criterion=env$criterion)
#           }
# )

# Seasonality auxiliary functions ----
# resolution: time steps in one cicle (e.g. 12 for months, 365 for days)
#' @export
seasonPars<- function(seasonRange){
  seasonMean<- apply(seasonRange, 1, function(x) (min(x) + max(x))/2)
  seasonAmplitude<- apply(seasonRange, 1, function (x) max(x) - min(x))
  
  return (data.frame(seasonMean, seasonAmplitude))
}

#' @export
seasonRange<- function(seasonMean, seasonAmplitude){
  seasonPars<- data.frame(seasonMean, seasonAmplitude)
  seasonRange<- with(seasonPars,  data.frame(min=seasonMean - seasonAmplitude/2, max=seasonMean + seasonAmplitude/2))
  return (seasonRange)
}

## Stochasticity
# beta distributded value

## Temporal correlation
# TODO

## Generic ----
#' @export
setMethod("print", signature(x="Env"),
          function(x, ...){
            print(S3Part(x), ...)
          }
)


#'  @export
setMethod("show", signature(object="Env"),
          function(object){
            cat("Object of class \"Environment\" with", nrow(object), "parameter sets\n")
            print(S3Part(object))
          }
)

# Only allowed to subset by rows but $ and [[i]] works for columns ## TODO needs a method Env(data.frame)
# `[.Env`<- function(x, ...){
#   Env(S3Part(x)[...])
# }

