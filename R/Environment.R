#' @include aaa-classes.R
NULL


## Constructor ----

#' @rdname Env
#'
#' @param pars 
#' @param seasonAmplitude 
#' @param seasonRange 
#' @param varJ 
#' @param varA
#' @param breedFail 
#'
#' @return a Env class object
#' @export
setGeneric("Env", function(pars, seasonAmplitude=c(0,1), seasonRange, varJ=c(0, 0.05), varA=c(0, 0.05), breedFail=c(0, 0.5, 0.9)) standardGeneric("Env"))

setMethod("Env",
          signature(pars="missing", seasonAmplitude="ANY", seasonRange="missing", varJ="ANY", varA="ANY", breedFail="ANY"),
          function(seasonAmplitude=c(0,1), varJ=c(0, 0.05), varA=c(0, 0.05), breedFail=c(0, 0.5, 0.9)){
            season<- getSeasonalParams(seasonAmplitude=seasonAmplitude, envMax=1)
            
            ## If seasonAmplitude == 0 there is no seasonality and the environment is constant (P_I(env) = P_i * env$seasonalMean)
            # season$seasonalMean[season$seasonAmplitude == 0]<- 1
            # season<- unique(season)
            
            comb<- expand.grid(varJ=varJ, varA=varA, breedFail=breedFail)
            comb<- merge(season, comb)
            
            comb<- unique(comb)
            comb<- data.frame(idEnv=formatC(seq_len(nrow(comb)), format="d", flag="0", width=nchar(nrow(comb))),
                              comb, stringsAsFactors=FALSE)
            
            new("Env", comb, seasonRange=getSeasonalRange(seasonMean=comb$seasonMean, seasonAmplitude=comb$seasonAmplitude))
          }
)

setMethod("Env",
          signature(pars="missing", seasonAmplitude="missing", seasonRange="matrix", varJ="ANY", varA="ANY", breedFail="ANY"),
          function(seasonRange, varJ=c(0, 0.05), varA=c(0, 0.05), breedFail=c(0, 0.5, 0.9)){
            mat<- matrix(seasonRange, ncol=2)
            seasonPars<- getSeasonalParams(seasonRange=mat)
            
            env<- Env(varJ=varJ, varA=varA, seasonAmplitude=seasonPars$seasonAmplitude, breedFail=breedFail)
            return (env)
          }
)  

setMethod("Env",
          signature(pars="data.frame", seasonAmplitude="missing", seasonRange="missing", varJ="missing", varA="missing", breedFail="missing"),
          function(pars){
            if (inherits(pars, "Model"))
              pars<- data.frame(pars, stringsAsFactors=FALSE)
            
            if (!"idEnv" %in% names(pars))
              pars$idEnv<- formatC(seq_len(nrow(pars)), format="d", flag="0", width=nchar(nrow(pars)))
            
            pars<- unique(pars[,c("idEnv", "seasonAmplitude", "seasonMean", "varJ", "varA", "breedFail")])
            
            # Sort rows by idEnv using naturalsort (slow) if installed and idEnv is non numeric
            suppressWarnings(idEnvnum<- as.numeric(pars$idEnv))
            if (anyNA(idEnvnum)){
              if (requireNamespace("naturalsort", quietly=TRUE)){
                pars<- pars[naturalsort::naturalorder(pars$idEnv),]
              }else{
                pars<- pars[order(pars$idEnv),]
              }
            }else{
              pars<- pars[order(idEnvnum),]
            }
            
            rownames(pars)<- pars$idEnv
            
            new("Env", pars, seasonRange=getSeasonalRange(seasonMean=pars$seasonMean, seasonAmplitude=pars$seasonAmplitude))
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
            timeSteps<- sin(timeSteps)
            
            pat<- t(apply(env[c("seasonAmplitude", "seasonMean")], 1, function(x){
              timeSteps * x["seasonAmplitude"] / 2 + x["seasonMean"]
            }))
            dimnames(pat)<- list(NULL, t=1:length(timeSteps))
            
            return (pat)
          }
)

setMethod("seasonalPattern", 
          signature(env="Env", resolution="ANY", cicles="ANY"),
          function(env, resolution=12, cicles=1){ 
            callGeneric(env=data.frame(env), resolution, cicles)
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
            callGeneric(env=data.frame(env), resolution, nSteps, interval, criterion)
          }
)


# Seasonality auxiliary functions ----
# resolution: time steps in one cicle (e.g. 12 for months, 365 for days)
#' @export
getSeasonalParams<- function(seasonAmplitude, seasonRange, envMax=1){
  if (!missing(seasonAmplitude)){
    seasonMean<- envMax - seasonAmplitude / 2
  }else if (!missing(seasonRange)){
    seasonMean<- apply(seasonRange, 1, function(x) (min(x) + max(x))/2)
    seasonAmplitude<- apply(seasonRange, 1, function (x) max(x) - min(x))
  }
  
  return (data.frame(seasonAmplitude, seasonMean))
}

#' @export
getSeasonalRange<- function(seasonMean, seasonAmplitude){
  seasonParams<- data.frame(seasonMean, seasonAmplitude)
  seasonRange<- with(seasonParams,  data.frame(min=seasonMean - seasonAmplitude / 2, max=seasonMean + seasonAmplitude / 2))
  
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


#' Env summary
#'
#' @rdname Env
#' @param object 
#' @param ... 
#'
#' @return
#' @export
summary.Env<- function(object, ...){
  lapply(S3Part(object)[, sapply(object, is.numeric)], function(x) sort(unique(x)))
}


#' Plot Env
#'
#' @rdname Env
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom graphics plot
plot.Env<- function(x, ...){
  x<- S3Part(x)

  cols<- intersect(names(x), c("seasonAmplitude", "seasonMean", "varJ", "varA", "breedFail"))
  x<- unique(x[, cols])
  out<- graphics::plot(x[, sapply(x, is.numeric)], ...)

  return(invisible(out))
}
