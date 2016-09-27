#' Find N0 for a given probability of persistence
#'
#' @name findN0_Pest
#' @param model a \code{Model} object.
#' @param Pobjective the probability that a given population still exists at the end of the simulations.
#' @param cl
#' @param verbose
#' @return for \code{findN0_Pest.scenario} a \code{data.frame} with N0, the probability to survive, 
#'   N0interpoled and objectiveprobability. For \code{findN0_Pest} a \code{Model} object containing the same 
#'   data.frame in the \code{model@sim@N0_Pest} slot.
#' @export
#'
#' @examples
findN0_Pest<- function(model=Model(), cl=parallel::detectCores(), Pobjective=.5, verbose=FALSE){
  scenario<- S3Part(model)
  scenario<- split(scenario, rownames(scenario))
  pars<- model@sim@params
  sim<- model@sim
  
  if (is.numeric(cl)){
    numericCL<- TRUE
    cl<- parallel::makeCluster(cl)
  } else {
    numericCL<- FALSE
  }
  
  parallel::clusterExport(cl=cl, c("sim", "Pobjective", "verbose"), envir=environment())
  parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
  parallel::clusterEvalQ(cl, library(LHR))

  N0_Pest<- parallel::parLapply(cl=cl, scenario, findN0_Pest.scenario, sim=sim, Pobjective=Pobjective, verbose=verbose)

  # N0_Pest<- lapply(scenario, findN0_Pest.scenario, sim=sim, Pobjective=Pobjective, verbose=verbose)
  # 
  # N0_Pest<- list()
  # for (i in seq_along(scenario)){
  #   N0_Pest[[i]]<- findN0_Pest.scenario(scenario=scenario[[i]], sim=sim, Pobjective=Pobjective, verbose=verbose)
  # }

  N0_Pest<- do.call("rbind", N0_Pest)
  
  simRes<- model@sim
  simRes@N0_Pest<- N0_Pest
  
  out<- Model(pars=S3Part(model), sim=simRes) # constructor removes results
  out@sim<- simRes

  if (numericCL) parallel::stopCluster(cl)
  
  return(out)
}


#' @rdname findN0_Pest
#' @param scenario a one row data.frame with the parameters.
#' @param sim a Sim object specifyng the simulation type and parameters.
#'
#' @return
#' @export
#'
#' @examples
findN0_Pest.scenario<- function(scenario=data.frame(Model(lh=LH(lambda=1))[1,]), sim=Sim.discretePopSim(replicates=10000),
                Pobjective=.5, verbose=FALSE){
  parsSim<- sim@params
  fun<- switch(class(sim),
              Sim.discretePopSim=Pestablishment.discretePopSim,
              Sim.numericDistri=Pestablishment.numericDistri,
              Sim.ABM=Pestablishment.ABM,
              Sim.ssa=Pestablishment.ssa)
  minN<- 1
  Pmin<- fun(N0=minN, scenario=scenario, parsSim=parsSim)
  
  if (inherits(sim, c("Sim.ABM", "Sim.ssa"))){
    if (!is.list(parsSim$N0)) parsSim$N0<- list(parsSim$N0)
    maxN<- max(sapply(parsSim$N0, sum))
    
    N0<- parsSim$N0[[1]]
    N0[N0 != 0]<- 1
    parsSim$N0<- list(N0)
    names(parsSim$N0)<- paste0("N", sum(N0))
    
    Pmax<- fun(N0=maxN, scenario=scenario, parsSim=parsSim)
  }else{
    maxN<- max(parsSim$N0)
    Pmax<- fun(N0=maxN, scenario=scenario, parsSim=parsSim)
  }
  
  N0<- N0anterior<- 0
  i<- 1
  found<- FALSE
  Pest<- NA
  
  while (!found){
    N0<- interpole(Pobjective, minN, maxN, Pmin, Pmax)
    N0<- ceiling(N0)
    
    if (is.na(N0)) break
      
    if (N0 < 1) {N0<- 1}
    
    if (N0 > parsSim$maxN){ # Evita passar Inf com a N0 (important per lambdes baixes)
      N0<- parsSim$maxN -1
      maxN<- parsSim$maxN
      Pmax<- 1
    } 
    
    if (N0 == N0anterior) {N0<- N0-1}
    
    # cat("N0=", N0, "\tNmin=", minN, "\tNp50max=", maxN, "\tPmin=", Pmin, "\tPmax=", Pmax, "\n")
    Pest<- fun(N0=N0, scenario=scenario, parsSim=parsSim)
    
    if (Pest < Pobjective){
      if (N0 < minN){
        maxN<- N0
        Pmax<- Pest
      }else{
        minN<- N0
        Pmin<- Pest
      } #Optimitzacio en cas que N0 estigui per sota de Pobjective, minN i maxN
    }else{# && maxN > N0){
      if (N0 > maxN){
        minN<-N0
        Pmin<- Pest
      }else{
        maxN<- N0
        Pmax<- Pest
      } #Optimitzacio en cas que N0 estigui per sobre de Pobjective, minN i maxN
    }
    if (maxN < minN){
      b<- maxN
      maxN<- minN
      minN<- b
      b<- Pmax
      Pmax<- Pmin
      Pmin<- b
    }
    
    if (verbose){cat("N0=", N0, "\tPestabliment=", Pest, "\tNmin=", minN, "\tNmax=", maxN, "\tPmin=", Pmin, "\tPmax=", Pmax, "\n")}
    
    if (maxN - minN < 2){
      if (maxN - minN == 1){
        if ((Pmax - Pobjective)^2 < (Pmin - Pobjective)^2){
          N0<- maxN
          Pest<- Pmax
        }else{
          N0<- minN
          Pest<- Pmin
        }
      }
      # cat("N0=", N0, "\tPestabliment=", Pest, "\n")
      found<- TRUE
    }
    N0anterior<- N0
    i<- i+1
  }
  N0interpoled<- interpole(Pobjective, minN, maxN, Pmin, Pmax)
  
  if (is.na(N0interpoled) & maxN == 1 & minN == 1){
      ## N0_Pest < 1
      N0interpoled<- interpole(Pobjective, 0, maxN, 0, Pmax)
  }
  
  return (data.frame(idScenario=scenario$idScenario, N0_Pest=N0, Pest, N0interpoled=N0interpoled, Pobjective, stringsAsFactors=FALSE))
}


interpole<- function(objective, x1, x2, y1, y2, type="lineal"){
  b<- (y2 - y1) / (x2 - x1) 
  a<- y1 - b * x1
  return ((objective - a) / b)
}


Pestablishment.discretePopSim<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  parsSim$N0<- round(N0)
  parsSim$raw<- FALSE
  parsSim$Ntf<- FALSE
  res<- runScenario.discretePopSim(scenario=scenario, pars=parsSim, verbose=FALSE)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

Pestablishment.numericDistri<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  parsSim$N0<- round(N0)
  parsSim$raw<- FALSE
  # parsSim$Ntf<- FALSE ## TODO Not set in constructor Sim.numericDistri
  res<- runScenario.numericDistri(scenario=scenario, pars=parsSim, verbose=FALSE)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

Pestablishment.ABM<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  selClass<- which(parsSim$N0[[1]] > 0)
  N0class<- parsSim$N0[[1]] * (N0 %/% sum(parsSim$N0[[1]])) # split N0 evenly between classes with N0 != 0
  
  mod<- N0 %% sum(parsSim$N0[[1]])
  if (mod > 0) N0class[selClass[1:mod]]<- N0class[selClass[1:mod]] + 1 # add the module to the firsts classes with N0 != 0
  
  parsSim$N0<- list(N0class)
  
  res<- runScenario.ABM(scenario=scenario, pars=parsSim)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

Pestablishment.ssa<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  N0class<- parsSim$N0[[1]] * N0 %/% sum(parsSim$N0[[1]]) # split N0 between classes with N0 != 0
  N0class[which(N0class > 0)[1]]<- N0class[which(N0class > 0)[1]] + N0 %% sum(N0class) # add the remaining to the first class with N0 != 0
  
  res<- exploreSSA(x0L=N0class, params=scenario, transitionMat=parsSim$transitionMat, rateFunc=parsSim$rateFunc, maxTf=parsSim$tf, replicates=parsSim$replicates,
                   raw=FALSE, discretePop=FALSE, finalPop=FALSE) # cl= TODO: pass as parameter
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}
