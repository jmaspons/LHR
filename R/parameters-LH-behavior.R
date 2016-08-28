# Parameters in probabilities. Check model-IBM-ssa_LH_behavior.R for ssa parameters.
# (get|set)(Params|Scenario|Behavior)

# Parameters ----

#' Title
#'
#' @param lh
#' @param env
#' @param strategies 
#' @param habDiffScenario
#' @param behavior 
#'
#' @return
#' @export
#'
#' @examples
getParamsCombination.LH_Beh<- function(lh=LH(), env=Env(seasonAmplitude=0, var=0), 
                                       habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"), 
                                       behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2")){
  habDiffScenario<- match.arg(habDiffScenario, several.ok=TRUE)
  behavior<- match.arg(behavior, several.ok=TRUE)
  
  if (any(env$seasonAmplitude > 0)){
    warning("Seasonality not implemented for this model. Scenarios discarded!")
    env<- env[env$seasonAmplitude == 0,]
  }
  
  lhEnv<- combineLH_Env(lh=lh, env=env)
  LH_Env<- lhEnv$scenario
  # parameters<- list(seasonBroodEnv=lhEnv$seasonBroodEnv) # TODO: add seasonality
  
  comb<- expand.grid(scenario=rownames(LH_Env), habDiffScenario=habDiffScenario, behavior=behavior, stringsAsFactors=FALSE)
  comb<- merge(comb, LH_Env, by.x="scenario", by.y=0)
  
  combL<- split(comb, rownames(comb))
  params<- lapply(combL, function (x){
    out<- getParams.LH_Beh(x, habDiffScenario=x$habDiffScenario, behavior=x$behavior)
    data.frame(scenario=x$scenario, habDiff=x$habDiffScenario, behavior=x$behavior, out, stringsAsFactors=FALSE)
  })
  params<- do.call(rbind, params)

  params<- merge(params, LH_Env[,setdiff(names(LH_Env), names(params))], by.x="scenario", by.y=0)
  rownames(params)<- with(params, paste(scenario, habDiff, behavior, sep="_"))
  
  # Sort
  params$scenario<- as.numeric(params$scenario)
  params$habDiff<- factor(params$habDiff, levels=c("identicalHab", "mortalHab2", "nestPredHab2"))
  params$behavior<- factor(params$behavior, levels=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2"))
    
  params<- params[order(params$scenario, params$habDiff, params$behavior),]
  
  return (params)
}


# returns a different strategies.and scenarios
### TODO: rename scenario -> lh. fix code
# diffHab2: named vector with the differences in the parameters at habitat 2 respect habitat 1
# Warning: clutch have no effect on the simulation. It's necessary to modify the reaction channels (getReactionChannels(clutch1, clutch2))
getParams.LH_Beh<- function(params=data.frame(b=1, broods=1, a=0.7, j=0.3, AFR=1, breedFail=0.5, jind=0.4615385, jbr=0.65),
                            diffHab2,
                            habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"),
                            behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2")){

  habDiffScenario<- match.arg(habDiffScenario)
  behavior<- match.arg(behavior)
  
  out<- with(params, data.frame(b1=b, b2=b,   broods=broods, PbF1=1 - jbr, PbF2=1 - jbr,  a1=a,ab1=a,j1=jind,  a2=a,ab2=a,j2=jind, AFR=AFR))
  
  # Add extra parameters for neutral behavior and density dependence
  out$K<- -1 # densodependence not implemented
  out$Pb1=1
  out$Pb2=1
  out$c1=1
  out$c2=1
  out$cF=1
  out$P1s=.5
  out$P1b=.5
  out$P1j=.5
  
  out<- split(out, rownames(out))
  
  if (!missing(diffHab2)){
    out<- lapply(out, function(x) setParams2diff1(x, diffHab2, type="probabilityMultiplicative"))
  }else{
    out<- lapply(out, function(x) setScenario(x, habDiffScenario, type="probabilityMultiplicative"))
    names(out)<- paste(names(out), habDiffScenario, sep="_")
    out<- lapply(out, function(x) setBehavior(x, behavior))
    names(out)<- paste(names(out), behavior, sep="_")
  }
  out<- data.frame(do.call("rbind", out))
  
  out<- lapply(out, function(x){
    attributes(x)<- NULL
    x
  })
  
  return (data.frame(out))
}


## TODO translate parameters to probabilities 
# Return new parames where parameters on habitat 2 are modified according to parameters from habitat 1 and a difference.
# diff: named vector with the proportion of difference respect habitat 1.
setParams2diff1<- function(params,
                           diff=getScenario("identicalHab"),
                           type=c("additive", "multiplicative", "probabilityMultiplicative", "lambda")){
  type<- match.arg(type)
  
  varsDiff<- gsub("Diff$", "", names(diff))
  misVars<- setdiff(c(paste0(varsDiff, 1), paste0(varsDiff, 2)), names(params))
  
  if (length(misVars) > 0){
    warning("Removing some variables from diff which are missing in params:\n", paste(misVars, collapse=", "))
    diff<- diff[!names(diff) %in% gsub("1", "Diff", misVars)]
    varsDiff<- gsub("Diff$", "", names(diff))
    if (length(diff) == 0) return (params)
  }
  
  selHab1<- sort(paste0(varsDiff, "1"))
  selHab2<- sort(paste0(varsDiff, "2"))
  
  diff<- diff[order(names(diff))]
  
  for (i in 1:length(diff)){
    params[selHab2[i]]<- switch(type,
                                additive= params[selHab1[i]] + diff[i],
                                multiplicative= params[selHab1[i]] * diff[i],
                                probabilityMultiplicative= 1 - (1 - params[selHab1[i]])^diff[i],
                                lambda="TODO")
  }
  
  return (params)
}

# params<- getParams.LH_Beh()
# setParams2diff1(params, diff=c(PbFDiff=.5, dDiff=-.3, gDiff=-.2))
getScenario<- function(habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2")){
  habDiffScenario<- match.arg(habDiffScenario)
  
  diff<- switch(habDiffScenario,
                `identicalHab`= c(bDiff=1, PbFDiff=1, aDiff=1, abDiff=1, jDiff=1),
                `mortalHab2`=   c(bDiff=1, PbFDiff=1, aDiff=0.5, abDiff=0.5, jDiff=0.5),
                `nestPredHab2`= c(bDiff=1, PbFDiff=2, aDiff=1, abDiff=1, jDiff=1)
  )
  return (diff)
}

setScenario<- function(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                       habDiffScenario="identicalHab", type="probabilityMultiplicative"){
  params<- setParams2diff1(params=params, diff=getScenario(habDiffScenario), type=type)
  
  return (params)
}


setBehavior<- function(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                       behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2")){
  behavior<- match.arg(behavior, several.ok=TRUE)
  
  if ("neutral" %in% behavior){
    params[c("Pb1","Pb2",  "c1", "c2", "cF",  "P1s","P1b","P1j")]<- c(1,1, 1,1,1, .5,.5,.5)
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
    params[c("P1s","P1b","P1j")]<- c(.8,.8,.8)
  }
  
  if ("preferHab2" %in% behavior){
    params[c("P1s","P1b","P1j")]<- c(.2,.2,.2)
  }

  
  return (params)
}


