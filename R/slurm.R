## run.slurm(): Simulate models using slurm cluster queue manager ----

#' @rdname Model
#'
#' @param model 
#' @param nodes
#' @param cpus_per_node  
#' @param ... parameters passed to \code{\link[rslurm]{slurm_apply}
#' @inheritParams rslurm::slurm_apply
#'
#' @return a \code{slurm_job} object to retrieve the results using \code{\link[rslurm]{get_slurm_out}}.
#' @examples 
#' sjob<- run.slurm(Model(), nodes=2, cpus_per_node=4, jobname="LHR", slurm_options=list(mem=100), submit=FALSE)
#' cleanup_files(sjob)
#' @include aaa-classes.R
#' @export
setGeneric("run.slurm", function(model, nodes, cpus_per_node, ...) standardGeneric("run.slurm"))

setMethod("run.slurm", 
          signature(model="Model", nodes="numeric", cpus_per_node="numeric"),
          function(model, nodes, cpus_per_node, ...){

            if (!requireNamespace("rslurm")) stop("run.slurm requires to install rslurm package.")
            
            nBatch<- nrow(model) %/% nodes
            
            splitScenarios<- suppressWarnings(split(1:nrow(model), rep(1:nBatch, times=nodes)))
            # sort(as.integer(unlist(splitScenarios)))
            
            params<- data.frame(i=seq_along(splitScenarios))
            
            call<- switch(class(model@sim),
                          Sim.discretePopSim="run.discretePopSim(model[splitScenarios[[i]],], cl=cpus_per_node)",
                          Sim.numericDistri="run.numericDistri(model[splitScenarios[[i]],], cl=cpus_per_node)",
                          Sim.ABM="run.ABM(model[splitScenarios[[i]],], cl=cpus_per_node, ...)")
            
            sjob<- rslurm::slurm_apply(f=s_call, params=params,
                                  add_objects=c("model", "splitScenarios", "call", "cpus_per_node"),
                                  nodes=nodes, cpus_per_node=cpus_per_node, ...)
            
            return(sjob)
          }
)



## findN0_Pest.slurm(): Simulate models using slurm cluster queue manager ----

#' Find N0 for a given probability of persistence
#'
#' @name findN0_Pest.slurm
#' @param model a \code{Model} object.
#' @param nodes
#' @param cpus_per_node
#' @param Pobjective the probability that a given population still exists at the end of the simulations.
#' @param verbose
#' @param ... parameters passed to \code{\link[rslurm]{slurm_apply}
#' @inheritParams rslurm::slurm_apply
#'
#' @return a \code{\link[rslurm]{slurm_job}} object to retrieve the results using \code{\link[rslurm]{get_slurm_out}}.
#' @examples 
#' sjob- findN0_Pest.slurm(Model(), nodes=2, cpus_per_node=4, jobname="LHR", slurm_options=list(mem=100), submit=FALSE)
#' res<- ensembleModel.slurm(sjob)
#' cleanup_files(sjob)
#' @include aaa-classes.R
#' @export
findN0_Pest.slurm<- function(model=Model(), nodes, cpus_per_node, Pobjective=.5, verbose=FALSE, ...){
  
  if (!requireNamespace("rslurm")) stop("run.slurm requires to install rslurm package.")
  
  nBatch<- nrow(model) %/% nodes
  
  splitScenarios<- suppressWarnings(split(1:nrow(model), rep(1:nBatch, times=nodes)))
  # sort(as.integer(unlist(splitScenarios)))

  params<- data.frame(i=seq_along(splitScenarios))
  
  call<- "findN0_Pest(model=model[splitScenarios[[i]],], cl=cpus_per_node, Pobjective=Pobjective, verbose=verbose)"
  
  sjob<- rslurm::slurm_apply(f=s_call, params=params,  
                     add_objects=c("model", "splitScenarios", "call", "cpus_per_node", "Pobjective", "verbose"),
                     nodes=nodes, cpus_per_node=cpus_per_node, ...)
  
  return(sjob)
}


## Function to get and ensamble the results ----
#' Recover results from simulations using slurm cluster
#'
#' @name ensembleModel.slurm
#' @param sjob a \code{slurm_job}
#' @inheritParams rslurm::get_slurm_out
#'
#' @return a \code{\link{Model}} object with all the results from a  \code{\link[rslurm]{get_slurm_out}}.
#' @examples 
#' sjob<- run.slurm(Model(), nodes=2, cpus_per_node=4,
#'          jobname="LHR", slurm_options=list(mem=100), submit=FALSE)
#' res<- ensembleModel.slurm(sjob)
#' cleanup_files(sjob)
#' @seealso \code{\link{run.slurm}}, \code{\link{findN0_Pest.slurm}}.
#' @include aaa-classes.R
#' @export
ensembleModel.slurm<- function(sjob, wait=TRUE){
  out<- rslurm::get_slurm_out(sjob, outtype="raw", wait=wait)
  # do.call(rbind, out)
  do.call(rbind, out)
}


## helpers ----

s_call<- function(i){
  eval(parse(text=call))
}

## TRUENO
# sopt<- list(partition="q-express", mem=100) # https://slurm.schedmd.com/sbatch.html

## After sending jobs
# sjob<- slurm_job(jobname="testRandom", nodes=3) # necessary to manually recreate one if the job was submitted in a different R session
# print_job_status(sjob)
# res<- get_slurm_out(sjob, outtype="raw")
