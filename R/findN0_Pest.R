#' Find N0 for a given probability of persistence
#'
#' @name findN0_Pest
#' @param model a \code{Model} object.
#' @param Pobjective the probability that a given population still exists at the end of the simulations.
#' @param cl
#' @param pb if \code{TRUE} and \link[pbapply]{pbapply} package is installed, show a progress bar. It 
#' increases the communication overhead between the main process and nodes / child processes. In most cases,
#' using \code{pbobtions(use_lb=TRUE)} improves the concurrency.
#' @param verbose if \code{TRUE}, print information of the search status for each iteration.
#' @param debug if \code{TRUE} run the simulations in a simple loop and print information about the state.
#' 
#' @return for \code{findN0_Pest.scenario} a \code{data.frame} with N0, the probability to survive, 
#'   N0interpoled and objectiveprobability. For \code{findN0_Pest} a \code{Model} object containing the same 
#'   data.frame in the \code{model@sim@N0_Pest} slot.
#' @export
#'
#' @examples
findN0_Pest<- function(model=Model(), cl=parallel::detectCores(), Pobjective=.5, pb=FALSE, verbose=FALSE, debug=FALSE){
  scenario<- S3Part(model)
  scenario<- split(scenario, rownames(scenario))
  pars<- model@sim@params
  sim<- model@sim
  
  if (!debug){
    if (is.numeric(cl)){
      if (.Platform$OS.type == "windows"){
        cl<- parallel::makePSOCKcluster(cl)
      }else{
        cl<- parallel::makeForkCluster(cl)
      }
      on.exit(parallel::stopCluster(cl))
    }
    
    parallel::clusterExport(cl=cl, c("sim", "Pobjective", "verbose"), envir=environment())
    parallel::clusterSetRNGStream(cl=cl, iseed=NULL)
    parallel::clusterEvalQ(cl, library(LHR))
  
    if (pb & requireNamespace("pbapply", quietly=TRUE)){
      N0_Pest<- pbapply::pblapply(scenario, findN0_Pest.scenario, sim=sim, Pobjective=Pobjective, verbose=verbose, cl=cl)
    }else{
      N0_Pest<- parallel::parLapply(cl=cl, scenario, findN0_Pest.scenario, sim=sim, Pobjective=Pobjective, verbose=verbose)
    }
    
  } else {
    N0_Pest<- list()
    for (i in seq_along(scenario)){
      message(i, "/", length(scenario))
      print(scenario[[i]], row.names=FALSE)
      eTime<- system.time(N0_Pest[[i]]<- findN0_Pest.scenario(scenario=scenario[[i]], sim=sim, Pobjective=Pobjective, verbose=verbose))
      message("Elapsed time: ", eTime["elapsed"])
    }
  }
  
  N0_Pest<- do.call("rbind", N0_Pest)
  
  simRes<- model@sim
  simRes@N0_Pest<- N0_Pest
  
  out<- Model(pars=S3Part(model), sim=simRes) # constructor removes results
  out@sim<- simRes

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
findN0_Pest.scenario<- function(scenario=data.frame(Model(lh=LH(lambda=1))[1,]),
                                sim=Sim.discretePopSim(replicates=10000), Pobjective=.5, verbose=FALSE){
  parsSim<- sim@params
  fun<- switch(class(sim),
              Sim.discretePopSim=Pestablishment.discretePopSim,
              Sim.numericDistri=Pestablishment.numericDistri,
              Sim.ABM=Pestablishment.ABM)

  
  if (inherits(sim, "Sim.ABM")){
    if (!is.list(parsSim$N0)) parsSim$N0<- list(parsSim$N0)
    selN0<- ceiling(length(parsSim$N0) / 2)
    maxN<- sum(parsSim$N0[[selN0]])
    
    N0<- parsSim$N0[[selN0]]
    N0[N0 != 0]<- 1
    parsSim$N0<- list(N0)
    names(parsSim$N0)<- paste0("N", sum(N0))
  }else{
    maxN<- parsSim$N0[ceiling(length(parsSim$N0) / 2)]
  }
  
  Pmax<- fun(N0=maxN, scenario=scenario, parsSim=parsSim)
  
  if (Pmax == 0){
    maxN<- parsSim$maxN
    Pmax<- fun(N0=maxN, scenario=scenario, parsSim=parsSim)
    
    if (Pmax < Pobjective){ # maxN == parsSim$maxN -> no solution
      warning("Pestablishment < Pobjective for maxN allowed in ", rownames(scenario))
      return (data.frame(idScenario=scenario$idScenario, N0_Pest=maxN, Pest=Pmax, N0interpoled=NA_real_, Pobjective, stringsAsFactors=FALSE))
    }
  }
  
  minN<- 1
  Pmin<- fun(N0=minN, scenario=scenario, parsSim=parsSim)
  
  if (Pmin > Pobjective){
    # warning("Pestablishment < Pobjective for N0=1 in ", rownames(scenario))
    maxN<- minN
    Pmax<- Pmin
    minN<- 0
    Pmin<- 0
    N0interpoled<- interpole(Pobjective, minN, maxN, Pmin, Pmax)
    
    return (data.frame(idScenario=scenario$idScenario, N0_Pest=maxN, Pest=Pmax, N0interpoled=N0interpoled, Pobjective, stringsAsFactors=FALSE))
  }
  
  N0<- N0anterior<- 0
  found<- FALSE
  Pest<- NA
  
  
  while (!found){
    N0<- interpole(Pobjective, minN, maxN, Pmin, Pmax)
    N0<- ceiling(N0)
    
    if (is.na(N0)) break
      
    if (N0 < 1) { N0<- 1 }
    
    if (N0 > parsSim$maxN){ # Evita passar Inf com a N0 (important per lambdes baixes)
      N0<- parsSim$maxN - 1
      maxN<- parsSim$maxN
      Pmax<- fun(N0=maxN, scenario=scenario, parsSim=parsSim)
      
      if (Pmax < Pobjective){ # maxN == parsSim$maxN -> no solution
        warning("Pestablishment < Pobjective for maxN allowed in ", rownames(scenario))
        return (data.frame(idScenario=scenario$idScenario, N0_Pest=maxN, Pest=Pmax, N0interpoled=NA_real_, Pobjective, stringsAsFactors=FALSE))
      }
    }

    if (N0 == N0anterior || N0 == maxN) {
      N0<- N0 - 1
    } else if (N0 == minN) N0<- N0 + 1
    
    Pest<- fun(N0=N0, scenario=scenario, parsSim=parsSim)
    
    if (Pest < Pobjective){
      if (N0 < minN){ 
        # Optimitzation when Pest < Pobjective and N0 < minN < maxN
        # Pobjective < Pmin < Pmax
        maxN<- N0
        Pmax<- Pest
      }else{ # minN < N0 < maxN
        minN<- N0
        Pmin<- Pest
      } 
    }else{
      if (N0 > maxN){ 
        # Optimitzation when Pest > Pobjective and N0 > maxN > minN
        # Pobjective > Pmax > Pmin 
        minN<-N0
        Pmin<- Pest
      }else{ # minN < N0 < maxN
        maxN<- N0
        Pmax<- Pest
      }
    }
    
    if (maxN < minN){
      tmp<- maxN
      maxN<- minN
      minN<- tmp
      
      tmp<- Pmax
      Pmax<- Pmin
      Pmin<- tmp
    }

    if (verbose)
      message("Scenario ", rownames(scenario), " N0=", N0, "\tPestabliment=", Pest, "\tNmin=", minN, "\tPmin=", Pmin, "\tNmax=", maxN, "\tPmax=", Pmax)
    
    if (Pmax < Pmin){
      warning("Inestability in Pest. Increase the number of replicates or check errors in scenario", rownames(scenario))
      found<- TRUE
    }
    
    
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
      
      found<- TRUE
    }
    N0anterior<- N0
  } # End find N0 loop
  
  N0interpoled<- interpole(Pobjective, minN, maxN, Pmin, Pmax)
  
  if (is.na(N0interpoled)){
    if (maxN == 1 & minN == 1){
      ## N0_Pest < 1
      N0interpoled<- interpole(Pobjective, 0, maxN, 0, Pmax)
    }else if (minN == maxN){
      N0interpoled<- minN
      Pest<- mean(Pmin, Pmax)
    }else if (Pmax == Pmin){
      Pest<- Pmin
      N0interpoled<- mean(maxN, minN)
      warning("Inestability in Pest. Increase the number of replicates")
    }else if (debug){
      browser()
    }
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
  res<- runScenario.discretePopSim(scenario=scenario, pars=parsSim)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

Pestablishment.numericDistri<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  parsSim$N0<- round(N0)
  parsSim$raw<- FALSE
  # parsSim$Ntf<- FALSE ## TODO Not set in constructor Sim.numericDistri
  res<- runScenario.numericDistri(scenario=scenario, pars=parsSim)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

Pestablishment.ABM<- function(N0, scenario, parsSim){
  if (N0 == 0) return(0)
  
  selClass<- which(parsSim$N0[[1]] > 0)
  N0class<- parsSim$N0[[1]] * (N0 %/% sum(parsSim$N0[[1]])) # split N0 evenly between classes with N0 != 0
  
  mod<- N0 %% sum(parsSim$N0[[1]])
  if (mod > 0){
    N0class[selClass[1:mod]]<- N0class[selClass[1:mod]] + 1 # add the module to the firsts classes with N0 != 0
    randomizeN0<- TRUE # avoid bias due to the order of the classes
  }else{
    randomizeN0<- FALSE
  }
  parsSim$N0<- list(N0class)
  
  res<- runScenario.ABM(scenario=scenario, pars=parsSim, randomizeN0=randomizeN0)
  
  Pest<- 1 - res$stats[,"extinct"]
  
  return(Pest)
}

