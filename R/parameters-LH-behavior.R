# Parameters in probabilities.
# (get|set)(Params|Scenario|Behavior)

# API Parameters----

#' getParamsCombination.LH_Beh
#'
#' Combine parameters of LH-Env scenarios with 2 patch habitat differences and behavior.
#' 
#' @param lh a \code{\link{LH}} object.
#' @param env a \code{\link{Env}} object.
#' @param habDiffScenario
#' @param behavior 
#' @param cl The number of cores to use or a cluster object (\code{\link[parallel]{makeCluster}} or 
#'   \code{\link[snow]{makeCluster} from \href{https://cran.r-project.org/web/packages/snow/index.html}{snow} package}) 
#' @param pb if \code{TRUE} and \link[pbapply]{pbapply} package is installed, show a progress bar.
#'
#' @return
#' @export
#'
#' @examples
getParamsCombination.LH_Beh<- function(lh=LH(), env=Env(seasonAmplitude=0, varJ=0, varA=0),
                                       habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"),
                                       behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2"),
                                       cl=parallel::detectCores(), pb=FALSE){
  habDiffScenario<- match.arg(habDiffScenario, several.ok=TRUE)
  behavior<- match.arg(behavior, several.ok=TRUE)
  
  if (any(env$seasonAmplitude > 0)){
    warning("Seasonality and environmental variation not implemented for this model. Scenarios discarded!")
    env<- env[env$seasonAmplitude == 0 & env$varA == 0 & env$varJ == 0, ]
  }
  
  if (is.numeric(cl)){
    if (.Platform$OS.type == "windows"){
      cl<- parallel::makePSOCKcluster(cl)
    }else{
      cl<- parallel::makeForkCluster(cl)
    }
    on.exit(parallel::stopCluster(cl))
  }
  
  
  lhEnv<- combineLH_Env(lh=lh, env=env)
  LH_Env<- lhEnv$scenario
  # parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv) # TODO: add seasonality
  
  comb<- expand.grid(idScenario=LH_Env$idScenario, habDiffScenario=habDiffScenario, behavior=behavior, stringsAsFactors=FALSE)
  comb<- merge(comb, LH_Env, by="idScenario")
  
  combL<- split(comb, rownames(comb))
  
  if (pb & requireNamespace("pbapply", quietly=TRUE)){
    params<- pbapply::pblapply(combL, function (x){
      out<- getParams.LH_Beh(x, habDiffScenario=x$habDiffScenario, behavior=x$behavior)
      out<- data.frame(idScenario=x$idScenario, idHabDiff=x$habDiffScenario, idBehavior=x$behavior, out, stringsAsFactors=FALSE)
      
      out
    }, cl=cl)
  }else{
    params<- parallel::parLapply(cl=cl, combL, function (x){
      out<- getParams.LH_Beh(x, habDiffScenario=x$habDiffScenario, behavior=x$behavior)
      out<- data.frame(idScenario=x$idScenario, idHabDiff=x$habDiffScenario, idBehavior=x$behavior, out, stringsAsFactors=FALSE)
      
      out
    })
  }

  params<- do.call(rbind, params)

  params<- merge(params, LH_Env[,c("idScenario", setdiff(names(LH_Env), names(params)))], by="idScenario")
  params$idScenario<- with(params, paste(idScenario, idHabDiff, idBehavior, sep="_"))
  rownames(params)<- params$idScenario
  
  # Sort
  params$idHabDiff<- factor(params$idHabDiff, levels=c("identicalHab", "mortalHab2", "nestPredHab2"))
  params$idBehavior<- factor(params$idBehavior, levels=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2"))
  params<- params[naturalsort::naturalorder(params$idScenario),]
  
  return (params)
}


#' getParamsCombination.LHEnv_2pactchBeh
#'
#' Combine parameters of LH-Env scenarios with 2 patch habitat differences and behavior.
#' 
#' @param lh a \code{\link{LH}} object.
#' @param env a \code{\link{Env}} object.
#' @param patchScenario a data.frame with the difference respect habitat 1 in different columns as in \code{\link{getPatchScenario}}.
#' @param nonBreedingSurv survival for non breeding adults as a proportion of adult survival (eg. 1.5 for a 50\% increase).
#' @param cl The number of cores to use or a cluster object (\code{\link[parallel]{makeCluster}} or 
#'   \code{\link[snow]{makeCluster} from \href{https://cran.r-project.org/web/packages/snow/index.html}{snow} package}) 
#' @param pb if \code{TRUE} and \link[pbapply]{pbapply} package is installed, show a progress bar.
#'
#' @return
#' @export
#'
#' @examples
getParamsCombination.LHEnv_2patchBeh<- function(lh=LH(), env=Env(seasonAmplitude=0, varJ=0, varA=0),
                                                patchScenario=getPatchScenario(habDiffIntensity=1.5, behaviorIntensity=2),
                                                nonBreedingSurv=1.5,
                                                cl=parallel::detectCores(), pb=FALSE){
  if (any(env$seasonAmplitude > 0)){
    warning("Seasonality and environmental variation not implemented for this model. Scenarios discarded!")
    env<- env[env$seasonAmplitude == 0 & env$varA == 0 & env$varJ == 0, ]
  }
  
  if (is.numeric(cl)){
    if (.Platform$OS.type == "windows"){
      cl<- parallel::makePSOCKcluster(cl)
    }else{
      cl<- parallel::makeForkCluster(cl)
    }
    on.exit(parallel::stopCluster(cl))
  }
  
  
  lhEnv<- combineLH_Env(lh=lh, env=env)
  LH_Env<- lhEnv$scenario
  
  aNonBreed<- 1 - (1 - LH_Env$a)^nonBreedingSurv
  LH_Env<- cbind(LH_Env, with(LH_Env, data.frame(b1=b, b2=b,  PbF1=1 - jbr, PbF2=1 - jbr,
                                a1=aNonBreed, ab1=a, sa1=s, j1=jind,
                                a2=aNonBreed, ab2=a, sa2=s, j2=jind)))
  
  # parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv) # TODO: add seasonality
  
  comb<- merge(LH_Env, patchScenario)
  comb$idScenario<- paste(comb$idScenario, comb$idHabDiff, comb$idBehavior, sep="_")
  
  diffCols<- grep("Diff$", names(comb), value=TRUE)
  diffCols<- setdiff(diffCols, "idHabDiff")
  otherCols<- setdiff(names(comb), diffCols)
  
  combL<- split(comb, 1:nrow(comb))
  
  if (pb & requireNamespace("pbapply", quietly=TRUE)){
    params<- pbapply::pblapply(combL, function (x){
      out<- getParams2diff1(params=x[, otherCols], diff=x[, diffCols], type="probabilityMultiplicative")
      return(out)
    }, cl=cl)
  }else{
    params<- parallel::parLapply(cl=cl, combL, function (x){
      out<- getParams2diff1(params=x[, otherCols], diff=x[, diffCols], type="probabilityMultiplicative")
      return(out)
    })
  }
  
  params<- do.call(rbind, params)
  
  # Sort columns
  ordCols<- c(names(lhEnv$scenario), "idHabDiff")
  ordCols<- c(ordCols, setdiff(names(params), ordCols))
  params<- params[, ordCols]
  
  # Sort rows
  params<- params[naturalsort::naturalorder(params$idScenario),]
  
  rownames(params)<- params$idScenario
  
  return (params)
}


## Helpers ----

# returns a different strategies.and scenarios
# diffHab2: named vector with the differences in the parameters at habitat 2 respect habitat 1
# @params nonBreedingSurv increase factor on survival for non breeding adults.
getParams.LH_Beh<- function(params=data.frame(b=1, broods=1, a=0.7, s=.6, j=0.3, AFR=1, breedFail=0.5, jind=0.4615385, jbr=0.65),
                            diffHab2, nonBreedingSurv=2,
                            habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"),
                            behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2")){

  habDiffScenario<- match.arg(habDiffScenario)
  behavior<- match.arg(behavior)
  
  aNonBreed<- 1 - (1 - params$a)^nonBreedingSurv
  
  out<- with(params, data.frame(b1=b, b2=b,   broods=broods, PbF1=1 - jbr, PbF2=1 - jbr,  a1=aNonBreed,ab1=a,sa1=s,j1=jind,  a2=aNonBreed,ab2=a,sa2=s,j2=jind, AFR=AFR))
  
  # Add extra parameters for neutral behavior and density dependence
  # out$K<- -1 # densodependence not implemented
  out$Pb1=1
  out$Pb2=1
  out$c1=1
  out$c2=1
  out$cF=1
  out$P1s=.5
  out$P1b=.5
  out$P1sa=.5
  out$P1j=.5
  
  out<- split(out, rownames(out))
  
  if (!missing(diffHab2)){
    out<- lapply(out, function(x) getParams2diff1(x, diffHab2, type="probabilityMultiplicative"))
  }else{
    out<- lapply(out, function(x) setHabScenario(x, habDiffScenario, type="probabilityMultiplicative"))
    names(out)<- paste(names(out), habDiffScenario, sep="_")
    out<- lapply(out, function(x) setBehScenario(x, behavior))
    names(out)<- paste(names(out), behavior, sep="_")
  }
  out<- data.frame(do.call("rbind", out))
  
  return (data.frame(out))
}


#' getParams2diff1
#'
#' @param params a data.frame that must 
#' @param diff a named vector with the difference respect habitat 1.
#' @param type a character vector with the type of difference applied
#'  \code{("additive", "multiplicative", "probabilityMultiplicative", "lambda")}.
#'
#' @details
#'  all names in \code{diff} must be in params. Trailing "Diff" in parameter names is removed
#' @return Return \code{params} where parameters on habitat 2 are modified according to parameters from habitat 1 and a difference.
#' @export
getParams2diff1<- function(params,
                           diff=getDiffHabScenario("mortalHab2"),
                           type=c("additive", "multiplicative", "probabilityMultiplicative", "lambda")){
  type<- match.arg(type)
  
  varsDiff<- gsub("Diff$", "", names(diff))
  varsDiffParams<- c(paste0(varsDiff, 1), paste0(varsDiff, 2))
  misVars<- setdiff(varsDiffParams, names(params))
  
  if (length(misVars) > 0){
    if (0 < length(selMisVars<- varsDiff[varsDiff %in% names(params)])){
      # add parameters with trailing 1 and 2 if missing (varsDiffParams)
      selMisVarsHabs<- paste0(rep(selMisVars,  each=2), rep(1:2, length(selMisVars)))
      params<- cbind(params, structure(params[, rep(selMisVars,  each=2)], names=selMisVarsHabs))
    }
    if (misVars<- setdiff(varsDiffParams, names(params)) > 0){
      warning("Removing some variables from diff which are missing in params:\n", paste(misVars, collapse=", "))
      diff<- diff[!names(diff) %in% gsub("1", "Diff", misVars)]
      varsDiff<- gsub("Diff$", "", names(diff))
      
      if (length(diff) == 0) return (params)
    }
  }
  
  selHab1<- sort(paste0(varsDiff, "1"))
  selHab2<- sort(paste0(varsDiff, "2"))
  
  diff<- diff[order(names(diff))]
  
  for (i in 1:length(diff)){
    params[, selHab2[i]]<- switch(type,
                                additive= params[, selHab1[i]] + diff[i],
                                multiplicative= params[, selHab1[i]] * diff[i],
                                probabilityMultiplicative= 1 - (1 - params[, selHab1[i]])^diff[i],
                                lambda="TODO")
  }
  
  return (params)
}

# PatchScenario = habitatDiff + behavior
getPatchScenario<- function(habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"), habDiffIntensity=2,
                                behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "static", "preferHab1", "preferHab2"),
                                behaviorIntensity=2){
  ## Combine habitatDiffScenario and intensity
  habComb<- expand.grid(habDiffScenario=habDiffScenario, habDiffIntensity=habDiffIntensity, stringsAsFactors=FALSE)
  
  # Remove combinations with no effect
  if (any(sel<- habComb$habDiffScenario == "identicalHab")){
    selRm<- which(sel)[-1]
    
    if (length(selRm) > 0){
      habComb<- habComb[-selRm,]
    }
  }
  
  
  habDiff<- lapply(split(habComb, 1:nrow(habComb)), function(x){
      getDiffHabScenario(habDiffScenario=x$habDiffScenario, intensity=x$habDiffIntensity)
    })
  habDiff<- as.data.frame(do.call(rbind, habDiff))
  rownames(habDiff)<- ifelse(habComb$habDiffScenario == "identicalHab", "identicalHab",
                             paste0(habComb$habDiffScenario, "X", habComb$habDiffIntensity))
  
  if (any(duplicated(habDiff))){
    warning("Some duplicated behavior parametres witch shouldn't")
    habDiff<- unique(habDiff)
  }
  
  ## Combine behavior and intensity
  behComb<- expand.grid(intensity=behaviorIntensity, behavior=behavior, stringsAsFactors=FALSE)
  
  # Remove combinations with no effects
  selComb<- !(behComb$behavior  %in% c("skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2") & behComb$intensity == 1)
  selComb<- selComb & !(behComb$behavior %in% c("neutral", "static") & duplicated(behComb$behavior))
  
  behComb<- behComb[selComb,]
  
  beh<- lapply(split(behComb, 1:nrow(behComb)), function(x){
    getBehavior(behavior=x$behavior, intensity=x$intensity)
  })
  beh<- as.data.frame(do.call(rbind, beh))
  rownames(beh)<- ifelse(behComb$behavior %in% c("neutral", "static"), behComb$behavior, paste0(behComb$behavior, "X", behComb$intensity))
  
  # beh[sort(rownames(beh)),]
  if (any(duplicated(beh))){
    warning("Some duplicated behavior parametres witch shouldn't")
    beh<- unique(beh)
  }
  
  habDiff<- cbind(idHabDiff=rownames(habDiff), habDiff, stringsAsFactors=FALSE)
  beh<- cbind(idBehavior=rownames(beh), beh, stringsAsFactors=FALSE)
  
  patchScenario<- merge(habDiff, beh)
  rownames(patchScenario)<- paste0(patchScenario$idHabDiff, "_", patchScenario$idBehavior)

  return(patchScenario)
}


getDiffHabScenario<- function(habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"), intensity=2){
  habDiffScenario<- match.arg(habDiffScenario)
  
  diff<- switch(habDiffScenario,
                `identicalHab`= c(bDiff=1, PbFDiff=1, aDiff=1, abDiff=1, saDiff=1, jDiff=1),
                `mortalHab2`  = c(bDiff=1, PbFDiff=1, aDiff=1 / intensity, abDiff=1 / intensity, saDiff=1 / intensity, jDiff=1 / intensity),
                `nestPredHab2`= c(bDiff=1, PbFDiff=intensity, aDiff=1, abDiff=1, saDiff=1, jDiff=1)
  )
  return (diff)
}


getBehavior<- function(behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "static", "preferHab1", "preferHab2"),
                       intensity=2){
  behavior<- match.arg(behavior)
  
  ## Modifiers of probabilities
  P.pos<- 1 - 1 / intensity
  P.neg<- 1 / intensity
  P.mult<- 1 - (1 - .5)^intensity ## Preferences. Neutral = .5
  
  # init with neutral
  # params[c("Pb1","Pb2",  "c1", "c2", "cF",  "P1s","P1b","P1sa","P1j")]<- c(1,1, 1,1,1, .5,.5,.5,.5)
  params<- data.frame(Pb1=1, Pb2=1,  c1=1, c2=1, cF=1, P1s=.5, P1b=.5,  P1sa=.5, P1j=.5, row.names=behavior)
  
  if ("skip" %in% behavior){ ## probability to skip breeding on habitat selection
    params["Pb2"]<- P.neg
  }
  
  if ("learnBreed" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0, 0, P.pos)
  }
  
  # Avoid habitat 2 after exploring or breeding fail (1 timestep memory only)
  if ("learnExploreBreed" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0, P.pos, P.pos)
  }
  
  if ("static" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0,0,0)
  }
  
  if ("preferHab1" %in% behavior){
    params[c("P1s","P1b","P1sa","P1j")]<- P.mult
  }
  
  if ("preferHab2" %in% behavior){
    params[c("P1s","P1b","P1sa","P1j")]<- 1 - P.mult
  }
  
  return (params)
}

## DEPRECATED ----
setHabScenario<- function(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.8,ab1=.7,sa1=.6,j1=.25,  a2=.8,ab2=.7,sa2=.6,j2=.25, AFR=1, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                       habDiffScenario="identicalHab", type="probabilityMultiplicative"){
  params<- getParams2diff1(params=params, diff=getDiffHabScenario(habDiffScenario), type=type)
  
  return (params)
}


setBehScenario<- function(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.8,ab1=.7,sa1=.6,j1=.25,  a2=.8,ab2=.7,sa2=.6,j2=.25, AFR=1, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                       behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "static", "preferHab1", "preferHab2")){
  behavior<- match.arg(behavior, several.ok=TRUE)
  
  if ("neutral" %in% behavior){
    params[c("Pb1","Pb2",  "c1", "c2", "cF",  "P1s","P1b","P1sa","P1j")]<- c(1,1, 1,1,1, .5,.5,.5,.5)
  }
  
  if ("skip" %in% behavior){ ## Increase breeding costs or better habitat selection
    params[c("Pb1", "Pb2")]<- c(1, .2)
  }
  
  if ("learnBreed" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0, 0, .8)
  }
  
  # Avoid habitat 2 after exploring or breeding fail (1 timestep memory only)
  if ("learnExploreBreed" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0, .8, .8)
  }
  
  if ("static" %in% behavior){
    params[c("c1", "c2", "cF")]<- c(0,0,0)
  }
  
  if ("preferHab1" %in% behavior){
    params[c("P1s","P1b","P1sa","P1j")]<- c(.8,.8,.8,.8)
  }
  
  if ("preferHab2" %in% behavior){
    params[c("P1s","P1b","P1sa","P1j")]<- c(.2,.2,.2,.2)
  }
  
  return (params)
}


