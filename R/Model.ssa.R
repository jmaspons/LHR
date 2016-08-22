#' @include aaa-classes.R
NULL

run.ssa.deterministic<- function(model, cl=parallel::detectCores(), ...){
  x0L<- model@sim@params$N0
  params<- S3Part(model)
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }

  res<- exploreSSA.deterministic(params=params, transitionMat=model@sim@params$transitionMat, rateFunc=model@sim@params$rateFunc, cl=cl, ...)
  ## TODO: return the Model adding this result to model@sim@deterministic
  res<- new("Sim.ssa", res$stats, params=model@sim@params, raw=unclass(res))
  
  if (numericCL) parallel::stopCluster(cl)
  
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
  ggplot2::ggplot(data=res, ggplot2::aes(N0, 1 - decrease, colour=behavior)) + ggplot2::geom_line() + ggplot2::geom_point() + 
    ggplot2::facet_grid(environment ~ LH) + ggplot2::coord_cartesian(ylim=c(0, 1))
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
  mRes<- reshape2::melt(mRes, value.name="N")
  mRes$quantile<- as.numeric(gsub("V", "", mRes$variable))
  mRes$quantile<- mRes$quantile / length(unique(mRes$quantile))
  resSel<- mRes[which(mRes$N0 == 2),]
  ggplot2::ggplot(data=resSel, ggplot2::aes(quantile, N, colour=behavior, group=scenario)) + ggplot2::geom_line() + 
    ggplot2::facet_grid(environment ~ LH) # + ggplot2::geom_point()# + ggplot2::coord_cartesian(ylim=c(0, 1))
}

##TODO: plot grow rates ~ time
#' @rdname Model.ssa
#'
#' @param gr ???
#' @param firstOnly 
#' @param omitOutliers 
#' @param ylim 
#'
#' @export
plotR<- function(gr, firstOnly=TRUE, omitOutliers=TRUE, ylim=stats::quantile(gr$r, probs=c(.025, .975), na.rm=TRUE, names=FALSE)){
  if (firstOnly) gr<- gr[grep("[", gr$period, fixed=TRUE),]
  gg<- ggplot2::ggplot(data=gr, ggplot2::aes(N0, r, colour=behavior)) + 
    ggplot2::facet_grid(environment ~ LH) + ggplot2::coord_cartesian(ylim=ylim) + ggplot2::geom_hline(yintercept=0, colour="red")
  if (omitOutliers){
    gg + ggplot2::geom_boxplot(outlier.size=0)
  }else{
    gg + ggplot2::geom_boxplot()
  }
}

#' @rdname Model.ssa
#'
#' @param gr ????
#' @param firstOnly 
#' @param omitOutliers 
#' @param ylim 
#'
#' @export
plotLambda<- function(gr, firstOnly=TRUE, omitOutliers=TRUE, ylim=stats::quantile(gr$lambda, probs=c(.025, .975), na.rm=TRUE, names=FALSE)){
  if (firstOnly) gr<- gr[grep("[", gr$period, fixed=TRUE),]
  gg<- ggplot2::ggplot(data=gr, ggplot2::aes(N0, lambda, colour=behavior)) + ggplot2::facet_grid(environment ~ LH) + 
    ggplot2::coord_cartesian(ylim=ylim) + ggplot2::geom_hline(yintercept=1, colour="red")
  if (omitOutliers){
    gg + ggplot2::geom_boxplot(outlier.size=0)
  }else{
    gg + ggplot2::geom_boxplot()
  }
}
