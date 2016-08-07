#' Model.ssa
#' 
#' @name Model.ssa
#'
#' @include Model.R
#' @export
setClass("Model.ssa", contains="Model")

## Model.ssa ----
#' @rdname Model.ssa
#'
#' @param pars 
#' @param sim 
#'
#' @return a \code{Model.ssa} object.
#' @examples Model.ssa()
#' 
#' @export
setGeneric("Model.ssa", function(pars=getParams.LH_Beh(), sim=Sim.ssa()) standardGeneric("Model.ssa"))

setMethod("Model.ssa",
          signature(pars="data.frame", sim="Sim.ssa"),
          function (pars, sim){
            new("Model.ssa", pars, sim=sim)
          }
)

setMethod("Model.ssa",
          signature(pars="ANY", sim="ANY"),
          function (pars=getParams.LH_Beh(), sim=Sim.ssa()){
            new("Model.ssa", pars, sim=sim)
          }
)

## Called from run(Model.ssa)
run.ssa<- function(model, cl=makeCluster(cores, type="FORK"), cores=detectCores(), ...){
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
  #   cores=1
  #   mc.preschedule=TRUE
  
  res<- exploreSSA(x0L=x0L, params=params, transitionMat=transitionMat, rateFunc=rateFunc, 
                   tf=tf, replicates=replicates, discretePop=discretePop, finalPop=finalPop, cl=cl, ...)
  res<- new("Sim.ssa", res$stats, Ntf=res$Ntf, params=model@sim@params, raw=unclass(res))
  return (res)
}

run.ssa.deterministic<- function(model, cl=makeCluster(cores, type="FORK"), cores=detectCores(), ...){
  x0L<- model@sim@params$N0
  params<- S3Part(model)

  res<- exploreSSA.deterministic(params=params, transitionMat=model@sim@params$transitionMat, rateFunc=model@sim@params$rateFunc, cl=cl, ...)
  ## TODO: return the Model adding this result to model@sim@deterministic
  res<- new("Sim.ssa", res$stats, params=model@sim@params, raw=unclass(res))
  return (res)
}


## Graphics ----
#' @rdname Model.ssa
#'
#' @param model 
#'
#' @export
plotPest<- function(model){
  res<- result(model)
  ggplot(data=res, aes(N0, 1 - decrease, colour=behavior)) + facet_grid(environment ~ LH) + geom_line() + geom_point() + coord_cartesian(ylim=c(0, 1))
}

#' @rdname Model.ssa
#'
#' @param model 
#'
#' @export
plotNtf<- function(model){
  res<- model@sim@Ntf
  scenario<- strsplit(rownames(res), "_")
  scenario<- do.call("rbind", scenario)
  colnames(scenario)<- c("LH", "environment", "behavior", "N0")
  scenario<- data.frame(scenario=rownames(res), scenario, stringsAsFactors=TRUE)
  scenario$N0<- ordered(as.numeric(gsub("N", "", scenario$N0)))
#   rownames(res)<- NULL
  rownames(scenario)<- rownames(res)
  mRes<- data.frame(scenario, res, stringsAsFactors=TRUE)
  mRes<- melt(mRes, value.name="N")
  mRes$quantile<- as.numeric(gsub("V", "", mRes$variable))
  mRes$quantile<- mRes$quantile / length(unique(mRes$quantile))
  resSel<- mRes[which(mRes$N0 == 2),]
  ggplot(data=resSel, aes(quantile, N, colour=behavior, group=scenario)) + facet_grid(environment ~ LH) + geom_line()# + geom_point()# + coord_cartesian(ylim=c(0, 1))
}

##TODO: plot grow rates ~ time
#' @rdname Model.ssa
#'
#' @param model 
#'
#' @export
plotR<- function(gr, firstOnly=TRUE, omitOutliers=TRUE, ylim=quantile(gr$r, probs=c(.025, .975), na.rm=TRUE, names=FALSE)){
  if (firstOnly) gr<- gr[grep("[", gr$period, fixed=TRUE),]
  gg<- ggplot(data=gr, aes(N0, r, colour=behavior)) + facet_grid(environment ~ LH) + coord_cartesian(ylim=ylim) + geom_hline(yintercept=0, colour="red")
  if (omitOutliers){
    gg + geom_boxplot(outlier.size=0)
  }else{
    gg + geom_boxplot()
  }
}

#' @rdname Model.ssa
#'
#' @param model 
#'
#' @export
plotLambda<- function(gr, firstOnly=TRUE, omitOutliers=TRUE, ylim=quantile(gr$lambda, probs=c(.025, .975), na.rm=TRUE, names=FALSE)){
  if (firstOnly) gr<- gr[grep("[", gr$period, fixed=TRUE),]
  gg<- ggplot(data=gr, aes(N0, lambda, colour=behavior)) + facet_grid(environment ~ LH) + coord_cartesian(ylim=ylim) + geom_hline(yintercept=1, colour="red")
  if (omitOutliers){
    gg + geom_boxplot(outlier.size=0)
  }else{
    gg + geom_boxplot()
  }
}
