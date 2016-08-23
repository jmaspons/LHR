#' @include aaa-classes.R
NULL


## Constructor ----

#' @rdname Env
#'
#' @param mean 
#' @param var 
#' @param seasonAmplitude 
#' @param seasonRange 
#' @param breedFail 
#'
#' @return a Env class object
#' @export
setGeneric("Env", function(pars, mean=c(1,.5), var=c(0, 0.1), seasonAmplitude=c(0,1), seasonRange, breedFail=c(0, 0.5, 1)) standardGeneric("Env"))

setMethod("Env",
          signature(pars="missing", mean="ANY", var="ANY", seasonAmplitude="ANY", seasonRange="missing", breedFail="ANY"),
          function(mean=c(1,.5), var=c(0, 0.1), seasonAmplitude=c(0,1), breedFail=c(0, 0.5, 1)){
            season<- data.frame(mean, seasonAmplitude)
            comb<- expand.grid(var=var, breedFail=breedFail)
            comb<- merge(season, comb)
            
            ## Remove combinations of mean and var out of the domain of the Beta distribution
            # betaPars<- fbeta(mean=comb$mean, var=comb$var)
            # 
            # if (any(!is.finite(betaPars$shape1) & comb$var != 0)){
            #     comb<- comb[is.finite(betaPars$shape1) | comb$var == 0,]
            #     warning("Some combinations of mean and var fall outside the parameter space of the beta distribution.")
            # }
             
            new("Env", comb, seasonRange=seasonRange(mean, seasonAmplitude))
          }
)


setMethod("Env",
          signature(pars="missing", mean="missing", var="ANY", seasonAmplitude="missing", seasonRange="matrix", breedFail="ANY"),
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
          signature(pars="data.frame", mean="ANY", var="ANY", seasonAmplitude="ANY", seasonRange="ANY", breedFail="ANY"),
          function(pars){
            pars<- unique(pars[,c("mean", "seasonAmplitude", "var", "breedFail")])
            
            new("Env", pars, seasonRange=seasonRange(pars$mean, pars$seasonAmplitude))
          }
)  

## Seasonal pattern ----
#' @rdname Env 
#'
#' @param env a \code{data.frame} of and \code{\link{Env}} object.
#' @param resolution 
#' @param cicles
#' 
#' @return A numeric vector with a sinusoidal pattern
#'
#' @export
setGeneric("seasonalPattern", function(env, resolution=12, cicles=1) standardGeneric("seasonalPattern"))


setMethod("seasonalPattern", 
          signature(env="data.frame", resolution="ANY", cicles="ANY"),
          function(env, resolution=12, cicles=1){ 
            timeSteps<- seq(0, 2*pi*cicles - 2*pi/resolution, by=2*pi/resolution)
            pat<- matrix(ncol=length(timeSteps), nrow=nrow(env), dimnames=list(NULL, t=1:length(timeSteps)))
            for (i in 1:nrow(env)){
              pat[i,]<- sin(timeSteps) * env$seasonAmplitude[i]/2 + env$mean[i]
            }
            return (pat)
          }
)

setMethod("seasonalPattern", 
          signature(env="Env", resolution="ANY", cicles="ANY"),
          function(env, resolution=12, cicles=1){ 
            callGeneric(env=S3Part(env), resolution, cicles)
          }
)

# Align a number of events separed by an interval to maximize a seasonal pattern of Env
#' @rdname Env
#'
#' @param env a \code{data.frame} or and \code{\link{Env}} object.
#' @param resolution the number of divisions of the sinusoidal pattern representing one year.
#' @param nSteps a vector with the number of events for each environent
#' @param interval 
#' @param criterion 
#'
#' @export
setGeneric("seasonOptimCal", function(env, resolution=12, nSteps=2, interval=1, criterion="maxFirst") standardGeneric("seasonOptimCal"))

setMethod("seasonOptimCal", 
          signature(env="data.frame", resolution="ANY", nSteps="ANY", interval="ANY", criterion="ANY"),
          function(env, resolution=12, nSteps=2, interval=1, criterion="maxFirst"){
            pattern<- seasonalPattern(env=env, resolution=resolution, cicles=1)
            n<- ncol(pattern) # resolution
            # Extend patterns 2 cicles to displace the event dates and find the best mean
            pattern<- matrix(rep(pattern,2), ncol=n * 2)
            
            # event calendar. SeasonAmplitude ensures that nrow(e) = nrow(env)
            e<- data.frame(nSteps, interval, season=env$seasonAmplitude)
            e$seasonLength<- e$nSteps * e$interval

            ## Maximize pattern on events
            if (criterion[1] == "maxFirst"){
              firstDate<- apply(pattern, 1, which.max)
              dates<- list()
              for (i in 1:nrow(e)){
                if (n < e$seasonLength[i]){
                  warning("Time steps don't fit in one cicle because number of steps (",
                          e$nSteps[i], ") and the interval (", e$interval[i], ") spans more than the resolution (", resolution, ") of one cicle.")
                  dates[[i]]<- NA
                }else{
                  dates[[i]]<- seq(firstDate[i], firstDate[i] + e$seasonLength[i], length=e$nSteps[i])
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

setMethod("seasonOptimCal", 
          signature(env="Env", resolution="ANY", nSteps="ANY", interval="ANY", criterion="ANY"),
          function(env, resolution=12, nSteps=2, interval=1, criterion="maxFirst"){
            callGeneric(env=S3Part(env), resolution, nSteps, interval, criterion)
          }
)


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
#' @rdname Env
#' @export
`[.Env`<- function(x, ...){
  Env(pars=data.frame(x)[...])
}

